import os, sys
import stat_collector as sc
import time
import argparse

GENERIC_SEQUENCE_TOKEN = '*'
OUTPUT_DELIMITER = '\t'

def count_unique_sequences(input_file, match_ranges, indexes=None):
    '''
    Returns a dictionary where each key corresponds to a set of sequences with
    the same nucleotides at the given ranges, and the value is a list of indexes
    in input_file where those sequences may be found. input_file may be any iterable
    of strings.
    If indexes is not None, it should be a set of indexes that the loop should
    check.
    '''

    input_file.seek(0)
    ret = {}
    for i, line in enumerate(input_file):
        if indexes is not None and i not in indexes:
            continue
        generic = get_generic_sequence(line.strip(), match_ranges)
        if generic in ret:
            ret[generic].add(i)
        else:
            ret[generic] = set([i])
    return ret

def get_generic_sequence(sequence, match_ranges):
    '''
    Returns a generic sequence where the bases outside match_ranges are denoted
    with the GENERIC_SEQUENCE_TOKEN symbol.
    '''
    match_ranges = sorted([(start % len(sequence), end % len(sequence)) for start, end in match_ranges])
    ret = ""
    for start, end in match_ranges:
        while len(ret) != 0 and len(ret) < start:
            ret += GENERIC_SEQUENCE_TOKEN
        ret += sequence[start:end]
    return ret

def write_hierarchical_unique_sequences(in_file, match_ranges, out_file, indexes=None, indent=0, complete_file=None, uniques=None):
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
        uniques = count_unique_sequences(in_file, [match_ranges[0]], indexes=indexes)
    for match, indexes in sorted(uniques.items(), reverse=True, key=lambda x: len(x[1])):
        if len(match_ranges) > 1:
            # Precompute the unique matches for the next match range to get the
            # number of unique items within this match
            new_uniques = count_unique_sequences(in_file, [match_ranges[1]], indexes=indexes)
            unique_submatches = len(new_uniques)
        else:
            new_uniques = None
            unique_submatches = len(indexes)

        string_to_write = "\t\t" * indent + OUTPUT_DELIMITER.join([str(len(indexes)), str(unique_submatches), match]) + "\n"

        if complete_file is not None:
            complete_file.write(string_to_write)
        if len(match_ranges) > 1:
            out_file.write(string_to_write)
            write_hierarchical_unique_sequences(in_file, match_ranges[1:], out_file, indexes=indexes, indent=indent + 1, complete_file=complete_file, uniques=new_uniques)

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

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    basename = os.path.basename(args.input)
    ranges = [(0, -15), (-15, -6)]
    with open(args.input, 'r') as file, open(os.path.join(args.output, basename), 'w') as out_file:
        if args.complete is not None:
            if not os.path.exists(args.complete):
                os.mkdir(args.complete)
            complete_file = open(os.path.join(args.complete, basename), 'w')
        else:
            complete_file = None
        write_hierarchical_unique_sequences(file, ranges, out_file, complete_file=complete_file)
        if complete_file is not None:
            complete_file.close()

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))
