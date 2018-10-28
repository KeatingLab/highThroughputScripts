## highThroughputScripts
**A pipeline for estimating Kds from cell surface display and high-throughput sequencing experiments.**

These scripts provide a flexible pipeline that begins with large high-throughput sequence read files and ends with a collection of protein sequences mapped to Kds. Each step of the process is designed to be modular and robust, with options to collect statistics about the runs. Furthermore, each step outputs the run parameters in a text file along with the actual output, making it easy to trace back variations between runs of the pipeline.

This repository contains two scripts from the TiteSeq paper, `KD_fit_log_poiss.py` and `ratio_x.py` (Adams et al 2016, [original on GitHub](https://github.com/jbkinney/16_titeseq)). Although the current work uses the Mass-Titer approach, these scripts can be executed in the `run_titeseq.py` step to provide an alternate method for fitting Kds.

Optionally, you can run the `random_sample_data.py` script to test the pipeline on a smaller version of the dataset. For example, to collect a random 10% of reads:

```
$ python3 random_sample_data.py forward_reads.fastq reverse_reads.fastq /output/dir/ -p 0.1
```

**Note: Other than `random_sample_data.py`, all scripts are written to be run in Python 2.**

### 1. Sort barcodes

The first and most time-intensive step, performed by `01_barcode_sort.py`, is simply to read the forward and reverse read files and sort the reads into the bins by their barcodes. Because this process is essentially limited by the machine's disk speed, this script may take up to about 30 hours to complete on a typical pair of Illumina fastq files (~60GB each). To run, simply call:

```
$ python 01_barcode_sort.py forward_reads.fastq reverse_reads.fastq /output/dir/
```

*Customization:* To adjust `01_barcode_sort.py` to fit a different barcoding scheme, the easiest approach is to change the `get_barcode_number` function or its associated constants. The current implementation computes a bin number based on the first five bases of the forward read and the first six bases of the reverse read, removing reads that do not meet a threshold of quality on these barcode regions.

### 2. Align reads

The script `02_align_reads.py` combines the forward and reverse reads by aligning both to an expected scaffold using a sliding-window approach. Type `python 02_align_reads.py -h` to learn about the various options available. This step is conducive to using a job manager such as SLURM:

```
#SBATCH --array=0-216
...

python 02_align_reads.py data/sorted_barcodes/barcode_${SLURM_ARRAY_TASK_ID} data/aligned_reads [-t 20] [-m 1] [-d 3] [-tr 20] [-mr 7] --stats [--check]
```

The `--check` parameter can be used to rerun the tasks partially, so that only the tasks for which output was not generated will be run. By default this script overwrites any files with the same names in the output directory.

*Customization:* You can adjust some parameters for the alignment at the top of the script. For example, the `REFERENCE_SEQUENCES`, `OUTPUT_RANGES`, and `REFERENCE_SCORING_RANGES` constants would be particularly useful to tailor to a particular experimental setup. (See the script for documentation on these constants.)

### 3. Count sequences

The script `03_count_sequences.py` counts unique sequences according to criteria you specify in the top of the script (see "Customization" below). It takes as input the contents of the `dnaframe` directory produced in step 2. To call this script using SLURM:

```
#SBATCH --array=0-216
...

input="data/aligned_reads/dnaframe/barcode_${SLURM_ARRAY_TASK_ID}_0"
output="data/sequence_counts"
python 03_count_sequences.py $input $output [-c $(output)_complete]
```

The result is a directory at `$output` with a tab-separated values file for each barcode, where each line contains a sequence and the number of counts for that sequence. If you specify multiple sorting tasks, the output also includes the number of unique sequences in the subsequent sorting tasks within that sequence.

*Customization:* The main parameters to adjust are `DISCARD_THRESHOLD` (number of mismatches to tolerate) and `SORTING_TASKS`, which indicates which bases the aligner should test. In most cases, only one sorting task will be needed. (The `-c` complete path command-line parameter is only necessary when you have multiple sorting tasks.)

### 4. Cluster and filter sequences

The primary function of `04_cluster_sequence_counts.py` is to further process the output of step 3 by condensing sequences that are off by a small number of mutations. Calling it is similar to step 3, and the format it returns is identical.

*Customization:* The `TASKS` variable indicates what jobs the script should perform, such as clustering by similarity, switching the hierarchical order for clustering (if multiple sorting tasks were used in step 3), and sorting the files by sequence count.

### 5. Summarize bin frequencies

At this point, we are ready to combine data from different bins to produce the distribution of each protein sequence across the fluorescence bins. The `05_write_titeseq_input.py` script tracks sequences across the set of bins you specify, and collects and normalizes them into a single CSV file. You can specify multiple tasks in the `TASKS` parameter to collect data from each experiment in parallel. Calling this script is as follows (here we process 4 experiments in parallel):

```
#SBATCH --array=0-3
...

input="data/collapsed_sequence_counts/task_1"
output="data/titeseq_input"
normpath="data/titeseq_sort_rates.csv"
python 05_write_titeseq_input.py $input $output -t ${SLURM_ARRAY_TASK_ID} -n $normpath
```

The purpose of the `$normpath` variable above is to normalize the counts on a bin-by-bin basis. This controls for having some bins oversampled and others undersampled during cell sorting. The `TASKS` parameter allows you to specify if you want to use this normalization procedure or not.

### 6. Fit dissociation curves

Finally, the `run_titeseq.py` script provides a few options for estimating Kds and fitting dissociation curves. It can be called simply by

```
#SBATCH --array=1-4
...

$ python run_titeseq.py data/titeseq_input_${SLURM_ARRAY_TASK_ID}.csv data/fit_output -m MODE
```

The `-m` (mode) parameter should be one of the following: `x_star`, `naive`, or `mass_titer`. The `x_star` mode uses the TiteSeq routine from Adams et al, 2016, while `naive` is a simple scikit-learn method and `mass_titer` uses lmfit. The output in all cases is a labeled CSV file ready for analysis.
