'''
The purpose of this script is to take a single file containing sequential pairs
of forward and reverse reads, and to align each pair together to produce a single
DNA sequence, protein sequence, and set of quality scores.
'''

from Bio import SeqIO,Seq
import multiprocessing
import sys, os
from itertools import izip, imap
import math
from functools import partial
import time
import argparse
from aligner import Aligner
import stat_collector as sc

FORMAT = 'fastq-sanger'
SEQUENCE_QUALITY_KEY = 'phred_quality'
UNSPECIFIED_BASE = 'N'

'''
These determine how statistics are collected when the --stats parameter is set.
The output will be in CSV format with three elements per row:

* The first number is a quality threshold (n / 10 gives the negative log probability
of a sequencing error).
* The second number is the number of bases that are allowed to fall below the
quality threshold.
* The third number is the number of reads that would pass the given threshold and
tolerance.
'''
STAT_FORWARD_KEY = "forward"
STAT_REVERSE_KEY = "reverse"
STAT_TOTAL_KEY = "total"
STAT_LENGTH_DIFF_KEY = "length_deltas"

STAT_DELETIONS_KEY = "deletions"
STAT_BAD_FORWARD_KEY = "bad forward read"
STAT_BAD_REVERSE_KEY = "bad reverse read"
STAT_EXCESS_LENGTH_KEY = "excessive length delta"

STAT_CUTOFFS = [0, 1, 2, 3, 5, 10, 20]
STAT_QUALITY_BINS = xrange(0, 40, 5)

'''
The references are the sequences onto which the forward and reverse reads will be
scaffolded. The output will be indexed by the number of the reference to which
each read was aligned.
'''
REFERENCE_SEQUENCES = ["ACTCGTCCCAAACAAGAACCTCAGGAAATCGATTTCCCGGACGATCTGCCAAACGACTTATCACG"]

'''
If OUTPUT_RANGES is not None, it should be a list with the same length as
REFERENCE_SEQUENCES. When a pair of reads is aligned to one of the reference
sequences, the sequence that is output will be trimmed to the range given by this
setting.
'''
OUTPUT_RANGES = None

'''
Same conditions as OUTPUT_RANGES, except a list of lists where each list corresponds
to one of the reference sequences. When a read is aligned to the reference sequence,
only the given ranges of the reference sequence will be used to score the alignment.
This is useful if leading or trailing regions are known to be varied.
'''
REFERENCE_SCORING_RANGES = None

'''
List of expression strings which will be evaluated and written out to the
params.txt file.
'''
PARAMETER_LIST = ["args.input", "args.output", "args.threshold", "args.misreads", "threshold_reverse", "misreads_reverse", "args.max_delta", "args.processes", "args.stats", "REFERENCE_SEQUENCES", "OUTPUT_RANGES", "REFERENCE_SCORING_RANGES", "STAT_CUTOFFS", "STAT_QUALITY_BINS"]

def combine_records(forward_record, reverse_record, reference_sequences, min_overlap=-1, max_overlap=-1, max_length_delta=1e30, reference_scoring_ranges=None):
    '''
    Computes the alignments of both forward and reverse reads to the reference
    sequences. Synthesizes those alignments, using the better-quality read in
    the case of a conflict. Returns (index, sequence, quality) where `index` is
    the index of the reference sequence used, `sequence` is the combined DNA
    sequence, and `quality` is the quality of each base in the combined sequence.

    The optional parameters min_overlap and max_overlap correspond to the overlap
    constraints on the alignment between the forward and reverse reads.
    '''
    aligner = Aligner()

    forward_str = str(forward_record.seq)
    reverse_str = str(reverse_record.seq.reverse_complement())

    # Align forward to references
    reference_index, forward_offset, forward_score = aligner.best_alignment(forward_str, reference_sequences, unidirectional=True, min_overlap=len(forward_str), candidate_scoring_ranges=reference_scoring_ranges)

    # Align forward to reverse
    reverse_offset, _ = aligner.align(forward_str, reverse_str, unidirectional=True, reverse=True, min_overlap=min_overlap, max_overlap=max_overlap)

    reference = reference_sequences[reference_index]
    reference_scoring_range = reference_scoring_ranges[reference_index] if reference_scoring_ranges is not None else None

    # Align reverse to reference
    reverse_offset_to_ref, reverse_score = aligner.align(reference, reverse_str, unidirectional=True, reverse=True, min_overlap=15, scoring_ranges=(reference_scoring_range, None))

    # Compare the pairwise scores of obeying the forward and obeying the reverse alignments to reference,
    # and adjust the alignment offsets accordingly.
    if reverse_score > forward_score:
        forward_offset = reverse_offset_to_ref - reverse_offset
        reverse_offset = reverse_offset_to_ref
    else:
        reverse_offset += forward_offset

    combined_sequence = ""
    combined_quality = []

    alignment_set = [(reference, 0), (forward_str, forward_offset), (reverse_str, reverse_offset)]
    # Uncomment to print the resulting alignments
    # print('\n'.join(aligner.format_multiple(*alignment_set)))

    # Discard the read if total length is too different from reference length
    if max_length_delta <= len(reference):
        if math.fabs(aligner.length(*alignment_set) - len(reference)) > max_length_delta:
            sc.counter(1, STAT_DELETIONS_KEY, STAT_EXCESS_LENGTH_KEY)
            return -1, None, None

    # Combine the reads to produce the overall sequence.
    # The aligner will enumerate the aligned characters or elements of each iterable we give it.
    # Zipping generators for both the sequence and the quality allows us to enumerate them together.
    sequence_generator = aligner.enumerate_multiple(*alignment_set)
    quality_generator = aligner.enumerate_multiple(([None for i in xrange(len(reference))], 0),
                                                   (forward_record.letter_annotations[SEQUENCE_QUALITY_KEY], forward_offset),
                                                   (reverse_record.letter_annotations[SEQUENCE_QUALITY_KEY], reverse_offset))
    for bases, qualities in izip(sequence_generator, quality_generator):
        _, forward_base, reverse_base = bases
        _, forward_quality, reverse_quality = qualities

        if forward_base is None and reverse_base is None:
            combined_sequence += UNSPECIFIED_BASE
            combined_quality.append(0)
        elif forward_base is None:
            combined_sequence += reverse_base
            combined_quality.append(reverse_quality)
        elif reverse_base is None:
            combined_sequence += forward_base
            combined_quality.append(forward_quality)
        else:
            base, quality = max([(forward_base, forward_quality), (reverse_base, reverse_quality)], key=lambda x: x[1])
            combined_sequence += base
            combined_quality.append(quality)

    return reference_index, combined_sequence, combined_quality

def combine_records_processor(references, (forward, reverse), threshold=0, misreads=0, threshold_reverse=0, misreads_reverse=0, output_ranges=None, **kwargs):
    '''
    Performs initial quality check and passes through to the main combine_records
    function. Returns (input, reference, dna_sequence, aa_sequence, quality), where reference
    is the index of the used reference. If an alignment was not computed because of
    bad quality, reference will be -1.

    If output_ranges is not None, its element corresponding to the chosen reference
    sequence will define the range of bases in the scaffolding sequence that
    should be output. (Note: it should match the reading frame of the coding
    sequence to ensure correct translation.)
    '''
    if not is_quality_read(forward, threshold=threshold, allowable_misreads=misreads):
        return (forward, reverse), -1, None, None, None, STAT_BAD_FORWARD_KEY
    if not is_quality_read(reverse, threshold=threshold_reverse, allowable_misreads=misreads_reverse):
        return (forward, reverse), -1, None, None, None, STAT_BAD_REVERSE_KEY

    reference, dna_sequence, quality = combine_records(forward, reverse, references, **kwargs)
    if reference == -1 or len(dna_sequence) == 0:
        return (forward, reverse), -1, None, None, None, STAT_EXCESS_LENGTH_KEY

    if output_ranges is not None:
        start, end = output_ranges[reference]
        dna_sequence = dna_sequence[start:end]
        quality = quality[start:end]

    aa_sequence = str(Seq.Seq(dna_sequence).translate())
    return (forward, reverse), reference, dna_sequence, aa_sequence, quality, None

### Helper functions

def is_quality_read(record, threshold=20, allowable_misreads=0):
    '''
    Checks the quality of the given SeqRecord. If more than `allowable_misreads`
    bases have quality below `threshold`, this function returns False, else it
    returns True.
    '''
    misread_count = sum(1 for q in record.letter_annotations[SEQUENCE_QUALITY_KEY] if q < threshold)
    return misread_count <= allowable_misreads

def update_quality_stats(parent_key, record):
    '''
    Helper function for write_combined_records that updates the given quality
    dictionary with qualities provided in the given SeqRecord (or list) object.
    '''
    qualities = record.letter_annotations[SEQUENCE_QUALITY_KEY] if type(record) == SeqIO.SeqRecord else record
    sorted_qualities = sorted(qualities)
    # Increment every (q, c) statistic for which the minimum quality
    # when the c worst bases are discarded is at least q.
    sc.permute_apply_counter([STAT_QUALITY_BINS, STAT_CUTOFFS], lambda x: 1 if sorted_qualities[x[1]] >= x[0] else 0, parent_key)

### Main function

def write_combined_records(input_path, references, out_dir, num_processes=15, stats=False, **kwargs):
    '''
    Combines the records at input_path by aligning them to the given reference sequences,
    and saves them to the appropriate locations within out_dir. The combined DNA
    sequences are written to out_dir/dnaframe, the amino acid sequences are
    written to out_dir/seqframe, and the qualities are written to qualframe (as CSV).

    The keyword arguments in **kwargs will be passed into the combine_records function.

    '''
    basename = os.path.basename(input_path)

    # Create directories if necessary
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    dna_path = os.path.join(out_dir, "dnaframe")
    aa_path = os.path.join(out_dir, "seqframe")
    qual_path = os.path.join(out_dir, "qualframe")
    for path in [dna_path, aa_path, qual_path]:
        if not os.path.exists(path):
            os.mkdir(path)
    if stats:
        stats_path = os.path.join(out_dir, "stats")

    # Open file streams
    dna_files = [open(os.path.join(dna_path, basename + "_{}".format(ref)), "w") for ref in xrange(len(references))]
    aa_files = [open(os.path.join(aa_path, basename + "_{}".format(ref)), "w") for ref in xrange(len(references))]
    qual_files = [open(os.path.join(qual_path, basename + "_{}".format(ref)), "w") for ref in xrange(len(references))]

    with open(input_path, 'rU') as file:

        records = SeqIO.parse(file, FORMAT)

        pool = multiprocessing.Pool(processes=num_processes)
        processor = partial(combine_records_processor,
                            references,
                            **kwargs)
        i = 0
        for result in pool.imap(processor, izip(records, records), chunksize=1000):
            if i % 1000 == 0:
                print("Finished {} sequences".format(i))
            i += 1
            original_input, ref_index, dna_sequence, aa_sequence, quality, error_key = result
            if stats:
                forward, reverse = original_input
                update_quality_stats(STAT_FORWARD_KEY, forward)
                update_quality_stats(STAT_REVERSE_KEY, reverse)
                if quality is not None:
                    update_quality_stats(STAT_TOTAL_KEY, quality)

                    delta = math.fabs(len(references[ref_index]) - len(dna_sequence))
                    sc.apply_counter(STAT_CUTOFFS, lambda c: delta <= c, STAT_LENGTH_DIFF_KEY)

            if ref_index == -1:
                if error_key is not None:
                    sc.counter(1, STAT_DELETIONS_KEY, error_key)
                continue

            dna_files[ref_index].write(dna_sequence + "\n")
            aa_files[ref_index].write(aa_sequence + "\n")
            qual_files[ref_index].write(",".join([str(q) for q in quality]) + "\n")

    for file in dna_files + aa_files + qual_files:
        file.close()
    if stats:
        sc.write(stats_path, prefix=basename)

if __name__ == '__main__':
    a = time.time()  # Time the script started
    parser = argparse.ArgumentParser(description='Reads a single barcode file containing pairs of forward and reverse reads, and aligns them to a set of possible scaffolds (defined at the top of this script). Writes the results to three directories, containing the DNA sequence, the AA sequence, and the quality scores, respectively.')
    parser.add_argument('input', metavar='I', type=str,
                        help='The path to the FASTQ input file')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the output directory')
    parser.add_argument('-t', '--threshold', type=int, default=0,
                        help='The quality score below which reads should not be used (default 0)')
    parser.add_argument('-m', '--misreads', type=int, default=0,
                        help='The number of allowable reads below the threshold score before the read is discarded')
    parser.add_argument('-tr', '--threshold-reverse', dest='threshold_reverse', type=int, default=-1,
                        help='The quality score below which reverse reads should not be used (default same as --threshold)')
    parser.add_argument('-mr', '--misreads-reverse', dest='misreads_reverse', type=int, default=-1,
                        help='The number of allowable reads below the threshold score before the reverse read is discarded (default same as --misreads)')
    parser.add_argument('-d', '--max-delta', dest='max_delta', type=int, default=1e30,
                        help='The maximum difference in length between the reference and the combined DNA alignment (default none)')
    parser.add_argument('-p', '--processes', type=int, default=15,
                        help='The number of processes to use')
    parser.add_argument('-s', '--stats', dest='stats', action='store_true',
                        help='Whether to output stats into output/[base]_stats.txt')
    parser.set_defaults(stats=False)
    args = parser.parse_args()

    threshold_reverse = args.threshold_reverse if args.threshold_reverse != -1 else args.threshold
    misreads_reverse = args.misreads_reverse if args.misreads_reverse != -1 else args.misreads

    write_combined_records(args.input,
                           REFERENCE_SEQUENCES,
                           args.output,
                           num_processes=args.processes,
                           threshold=args.threshold,
                           misreads=args.misreads,
                           threshold_reverse=threshold_reverse,
                           misreads_reverse=misreads_reverse,
                           min_overlap=10,
                           max_overlap=30,
                           stats=args.stats,
                           output_ranges=OUTPUT_RANGES,
                           max_length_delta=args.max_delta,
                           reference_scoring_ranges=REFERENCE_SCORING_RANGES)

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))

    # Write out the input parameters to file
    params_path = os.path.join(args.output, "params.txt")
    sc.write_input_parameters([(name, eval(name)) for name in PARAMETER_LIST], params_path)
