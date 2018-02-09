'''
This script takes the output of the 02_align_reads script that consists of lines
containing a single DNA sequence each, and produces a "sequence hierarchy" file
that expresses the counts of various fragments of sequence. For instance, if this
script is asked (in SORTING_TASKS) to count first by bases 0-10, then by bases
11-15, the output would be a text file with the following format:

5   3   AAAAAAAAAA
        2   2   BBBBB
        2   2   CCCCC
        1   1   DDDDD
4   1   EEEEEEEEEE
        4   4   FFFFF
...

The first number in each top-level item indicates the total number of reads that
had the given sequence in bases 0-10; the second number indicates the number of
unique sequences found in the range 11-15 that had this 0-10 sequence. A similar
pattern holds for the next-level items. This file format is used by the 04 and 05
scripts as well, so the tools can be swapped out as necessary. Tools to read and
write this file format can be found in seq_hierarchy_tools.py.
'''

import os, sys
import stat_collector as sc
import time
import argparse
from aligner import *

OUTPUT_DELIMITER = '\t'
STAT_SCORES = "scores"

'''
If the score of the alignment is less than the length of the template minus
this amount, the sequence will be discarded.
'''
DISCARD_THRESHOLD = 2

'''
Provides the criteria by which the input file will be sorted. If the task is a
tuple of integers (start, end), the script will search for the unique sequences
at those positions in each sequence.

If the task is a tuple of strings (template, output_template), the script will
attempt to align the template to each sequence in the file. The positions that
are denoted by VARIABLE_REGION_TOKEN in the output_template, which must be the
same length as template, will be assembled as the unique string for that sequence.
For example, if the sequence were ABCDEFGH, the template were AB**E*, and the
output template were A*****, the returned generic sequence would be '.BCDEF'.
'''
SORTING_TASKS = [(VARIABLE_REGION_TOKEN * 27 + "CCGGACGATCTGCC", VARIABLE_REGION_TOKEN * 41),
                 (-15, -6)]

'''
List of expression strings which will be evaluated and written out to the
params.txt file.
'''
PARAMETER_LIST = ["args.input", "args.output", "args.complete", "SORTING_TASKS", "DISCARD_THRESHOLD"]

def count_unique_sequences(input_file, match_task, indexes=None):
    '''
    If indexes is None, returns two items. First, a dictionary where each key
    corresponds to a set of sequences with the same nucleotides at the given
    position, and the value is a list of indexes in input_file where those
    sequences may be found. Second, the length of the file in lines.
    input_file may be any iterable of strings.

    If indexes is not None, it should be a list of the same length as the number
    of lines in the input. The return value in this case, will be a dictionary where
    each key is one of the values of the indexes list. The value at that key will
    be the unique sequences for all lines that have that given index value.
    '''

    input_file.seek(0)
    ret = {}
    num_lines = 0
    for i, line in enumerate(input_file):
        num_lines = i

        if indexes is not None:
            if indexes[i] is None:
                continue
            if indexes[i] not in ret:
                ret[indexes[i]] = {}
            ret_to_write = ret[indexes[i]]
        else:
            ret_to_write = ret

        if type(match_task[0]) is str:
            # Align the match_task template
            generic = get_generic_sequence_by_alignment(line.strip(), match_task[0], match_task[1])
        else:
            # Use the bases at the given positions
            generic = get_generic_sequence_by_position(line.strip(), [match_task])
        if generic is None:
            continue

        if generic in ret_to_write:
            ret_to_write[generic].add(i)
        else:
            ret_to_write[generic] = set([i])

    return ret, num_lines + 1

def get_generic_sequence_by_position(sequence, match_ranges):
    '''
    Returns a generic sequence where the bases outside match_ranges are denoted
    with the GENERIC_SEQUENCE_TOKEN symbol.
    '''
    match_ranges = sorted([(start % len(sequence), end % len(sequence)) for start, end in match_ranges])
    ret = ""
    ret_start = 0
    for start, end in match_ranges:
        if len(ret) == 0:
            ret_start = start
        else:
            while len(ret) < start - ret_start:
                ret += GENERIC_SEQUENCE_TOKEN
        ret += sequence[start:end]
    return ret

def get_generic_sequence_by_alignment(sequence, template, output_template):
    '''
    Returns a generic sequence by aligning template to sequence, and including
    only the bases in positions marked with VARIABLE_REGION_TOKEN in the template.

    For example, if the sequence were ABCDEFGH and the template were AB**E*, the
    returned generic sequence would be 'CD.F'.
    '''
    aligner = Aligner(different_score=0)
    offset, score = aligner.align(sequence, template, min_overlap=len(template))
    num_matching_bases = len([c for c in template if c not in NON_SCORED_TOKENS])
    sc.counter(1, STAT_SCORES, score)
    if score < num_matching_bases - DISCARD_THRESHOLD:
        return None

    ret = ""
    ret_start = 0
    for i, (base_1, base_2) in enumerate(aligner.enumerate(sequence, output_template, offset)):
        if base_2 == VARIABLE_REGION_TOKEN:
            if len(ret) == 0:
                ret_start = i
            else:
                while len(ret) < i - ret_start:
                    ret += GENERIC_SEQUENCE_TOKEN
            ret += base_1
    return ret

def write_hierarchical_unique_sequences(in_file, match_ranges, out_file, indent=0, complete_file=None, uniques=None, file_length=None):
    '''
    Groups the file in_file by the first match range, then each group by the
    second match range, and so on. Writes the distribution of each match to
    out_file (a file object).

    By default, the last group's keys will not be written to file. If desired,
    specify complete_file to write those keys in another hierarchical level to a
    separate file.
    '''
    if len(match_ranges) == 0:
        return
    if uniques is None:
        uniques, file_length = count_unique_sequences(in_file, match_ranges[0])

    sorted_uniques = sorted(uniques.items(), reverse=True, key=lambda x: len(x[1]))

    if len(match_ranges) > 1:
        # Precompute the unique matches for the next match range to get the
        # number of unique items within each match at this round
        index_list = [None for _ in xrange(file_length)]
        for i, (match, indexes) in enumerate(sorted_uniques):
            for index in indexes:
                index_list[index] = i
        new_uniques, _ = count_unique_sequences(in_file, match_ranges[1], indexes=index_list)
    else:
        new_uniques = None

    for i, (match, indexes) in enumerate(sorted_uniques):
        if len(match_ranges) > 1:
            unique_submatches = len(new_uniques[i])
        else:
            unique_submatches = len(indexes)

        string_to_write = "\t\t" * indent + OUTPUT_DELIMITER.join([str(len(indexes)), str(unique_submatches), match]) + "\n"

        if complete_file is not None:
            complete_file.write(string_to_write)
        if len(match_ranges) > 1:
            out_file.write(string_to_write)
            write_hierarchical_unique_sequences(in_file, match_ranges[1:], out_file, indent=indent + 1, complete_file=complete_file, uniques=new_uniques[i], file_length=file_length)

### Main function

def main_count_sequences(input, output, tasks, complete_path=None):
    if not os.path.exists(output):
        os.mkdir(output)
    basename = os.path.basename(input)
    with open(input, 'r') as file, open(os.path.join(output, basename), 'w') as out_file:
        if complete_path is not None:
            if not os.path.exists(complete_path):
                os.mkdir(complete_path)
            complete_file = open(os.path.join(complete_path, basename), 'w')
        else:
            complete_file = None

        write_hierarchical_unique_sequences(file, tasks, out_file, complete_file=complete_file)

        if complete_file is not None:
            complete_file.close()

    sc.write(os.path.join(output, "stats"), prefix=basename)

if __name__ == '__main__':
    a = time.time()  # Time the script started
    parser = argparse.ArgumentParser(description='Groups the sequences in the given file into hierarchies based on identical sequences in certain ranges.')
    parser.add_argument('input', metavar='I', type=str,
                        help='The path to the input file, where each line is a sequence')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the output directory')
    parser.add_argument('-c', '--complete', type=str, default=None,
                        help='The path to an additional output directory for the complete set of unique sequences')
    args = parser.parse_args()

    main_count_sequences(args.input, args.output, SORTING_TASKS, complete_path=args.complete)

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))

    # Write out the input parameters to file
    params_path = os.path.join(args.output, "params.txt")
    sc.write_input_parameters([(name, eval(name)) for name in PARAMETER_LIST], params_path)
