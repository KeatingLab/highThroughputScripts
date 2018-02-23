'''
This script takes a series of files containing collapsed sequence counts, tracks
the unique sequences through those files, and produces a CSV file that delineates
each sequence and its counts in each of the bins.
'''

from seq_hierarchy_tools import *
import time
import argparse
import stat_collector as sc
import os

TASK_TRACK_SEQS = "track_seqs"
TASK_TRACK_SINGLE_SEQ = "track_single_seq"

def generate_barcode_list(start, end):
    '''
    Convenience function that produces a list of filenames given barcode numbers.
    '''
    return ["barcode_{}_0".format(i) for i in range(start, end + 1)]

normalization_path = "../../titeseq_sort_rates.csv"
'''
If this parameter is nonzero, then the second half of the CSV numbers will be
enforced to be nonzero by adding this value to every element before normalization.
'''
NORMALIZATION_MINIMUM_VALUE = 0.1

'''
The tasks that the script should perform. The parameters in each tuple should be
as follows:
* If the task is TASK_TRACK_SEQS:
    * A list of file names that can be found in the input directory and should be
      tracked
    * The name of the output file to write to within the output directory
    * A path to a CSV file relative to the input directory containing
      normalization factors, by which each read count will be multiplied. The file
      should contain a line for each barcode file, where the first component is
      the filename and the second component is the multiplying factor. If this is
      None, no normalization will be applied.
* If the task is TASK_TRACK_SINGLE_SEQ:
    * A list of file names that can be found in the input directory and should be
      tracked
    * The individual sequence that should be searched for within each file
    * The name of the output file to write to within the output directory
'''
TASKS = [
    (TASK_TRACK_SEQS, generate_barcode_list(24, 55), "titeseq_input_1.csv", normalization_path),
    (TASK_TRACK_SEQS, generate_barcode_list(73, 104), "titeseq_input_2.csv", normalization_path),
    (TASK_TRACK_SEQS, generate_barcode_list(120, 151), "titeseq_input_3.csv", normalization_path),
    (TASK_TRACK_SEQS, generate_barcode_list(168, 199), "titeseq_input_4.csv", normalization_path)
]

'''
List of expression strings which will be evaluated and written out to the
params.txt file.
'''
PARAMETER_LIST = ["args.input", "args.output", "TASKS"]

def process_sequence_file(in_path, processed_counts, target_sequence=None, return_counts=False):
    '''
    in_path is a path to a sequence counts file, and processed_counts is a
    dictionary of sequences mapped to dictionaries of file names to counts. For
    instance:
    {'ABCDE': {'sequence_file_1': 3, 'sequence_file_2': 5, ...},
     'BCDEF': {'sequence_file_1': 1, 'sequence_file_2': 7, ...}}
    Returns the updated processed_counts dictionary.

    If return_counts is True, two additional items are returned: the total
    number of sequences in the file and the total number of reads in the file.
    '''
    basename = os.path.basename(in_path)

    with open(in_path, 'r') as file:
        seqs = read_sequence_dicts(file)

    total_reads = 0
    for seq in seqs:
        root, value = get_root_item(seq)
        count, unique = get_sequence_counts(value)
        total_reads += count

        if target_sequence is not None and root != target_sequence:
            continue

        if root not in processed_counts:
            processed_counts[root] = {}
        processed_counts[root][basename] = count

    if return_counts:
        return processed_counts, len(seqs), total_reads
    return processed_counts

def read_normalization_factors(path):
    '''
    Reads the normalization factor CSV file into a dictionary keyed by file name.
    '''
    ret = {}
    with open(path, 'r') as file:
        for line in file:
            comps = line.split(',')
            if len(comps) < 2: continue
            try:
                ret[comps[0]] = float(comps[1])
            except:
                continue
    return ret

def track_sequences(input_dir, output_dir, task):
    '''
    Tracks the unique sequences through the files given by the file names in
    `task`, and writes their counts in CSV format to the specified output file
    in the output directory.
    '''
    paths = task[1]
    out_path = os.path.join(output_dir, task[2])
    result = {}
    for path in paths:
        result = process_sequence_file(os.path.join(input_dir, path), result)
    print("Found {} unique sequences.".format(len(result)))

    normalization_path = os.path.join(input_dir, task[3])
    if normalization_path is not None:
        normalization_factors = read_normalization_factors(normalization_path)
    else:
        normalization_factors = {}

    with open(out_path, "w") as file:
        for seq, counts in result.items():
            comps_list = [seq]
            # First write all paths with no normalization
            for path in paths:
                if path in counts:
                    comps_list.append(str(counts[path]))
                else:
                    comps_list.append("0")
            # Now perform normalization if provided
            for path in paths:
                factor = normalization_factors[path] if path in normalization_factors else 1.0
                if path in counts:
                    comps_list.append(str((counts[path] + NORMALIZATION_MINIMUM_VALUE) * factor))
                else:
                    comps_list.append(str(NORMALIZATION_MINIMUM_VALUE * factor))
            file.write(",".join(comps_list) + "\n")

def track_single_sequence(input_dir, output_dir, task):
    '''
    Tracks a single sequence through the files given by the file names in
    `task`, and writes their counts in CSV format to the specified output file
    in the output directory.
    '''
    paths = task[1]
    target = task[2]
    out_path = os.path.join(output_dir, task[3])
    result = {}
    seq_counts = {}
    read_counts = {}
    for path in paths:
        result, seq_count, read_count = process_sequence_file(os.path.join(input_dir, path), result, target_sequence=target, return_counts=True)
        seq_counts[path] = seq_count
        read_counts[path] = read_count

    if len(result) == 0:
        print("Could not find sequence {} in {} sequence count files.".format(target, len(paths)))
        return
    seq, counts = result.items()[0]
    print("Found sequence {} in {} files.".format(target, len(counts)))
    with open(out_path, "w") as file:
        file.write(seq + "\n")
        file.write(",".join(["File Name", "Sequence Reads", "File Total Unique Sequences", "File Total Reads"]) + "\n")
        for path in paths:
            comps_list = [path]
            if path in counts:
                comps_list.append(str(counts[path]))
            else:
                comps_list.append("0")
            comps_list.append(str(seq_counts[path]))
            comps_list.append(str(read_counts[path]))
            file.write(",".join(comps_list) + "\n")

def perform_task(task, args):
    '''
    Delegates the given task to one of the functions in the script.
    '''
    if task[0] == TASK_TRACK_SEQS:
        track_sequences(args.input, args.output, task)
    elif task[0] == TASK_TRACK_SINGLE_SEQ:
        track_single_sequence(args.input, args.output, task)
    else:
        print("Unrecognized task type {}".format(task[0]))

if __name__ == '__main__':
    a = time.time()  # Time the script started
    parser = argparse.ArgumentParser(description='Tracks the sequences across sets of barcode hierarchy files according to the tasks listed at the top of the script. One file will be written out for each task.')
    parser.add_argument('input', metavar='I', type=str,
                        help='The path to the input sequence hierarchy directory')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the output directory')
    parser.add_argument('-t', '--task', type=int, default=-1,
                        help='The task index to perform (if using multiple nodes)')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    if args.task == -1:
        for i, task in enumerate(TASKS):
            print("Task {}...".format(i))
            perform_task(task, args)
    else:
        task = TASKS[args.task]
        perform_task(task, args)

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))

    # Write out the input parameters to file
    params_path = os.path.join(args.output, "params.txt")
    sc.write_input_parameters([(name, eval(name)) for name in PARAMETER_LIST], params_path)
