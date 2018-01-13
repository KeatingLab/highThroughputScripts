'''
This script simply concatenates the results of the statistics written by 02_align_reads.py.
Each file is in CSV format - see 02_align_reads.py for how the file is structured.
Note that this file will concatenate any set of distributions in CSV format, not just those
for 02_align_reads.py.
'''

import os
import argparse

SEPARATOR = '\t'

def combine_quality_stats(input_dir, marker, out_dir):
    quality_information = {}
    for path in os.listdir(input_dir):
        if path[0] == '.': continue
        if marker not in path: continue
        with open(os.path.join(input_dir, path), 'r') as file:
            for line in file:
                comps = line.strip().split(',')
                if len(comps) < 2:
                    print("Not enough elements to aggregate statistics.")
                    break
                current_dict = quality_information
                for comp in comps[:-2]:
                    if comp not in current_dict:
                        current_dict[comp] = {}
                    current_dict = current_dict[comp]
                if comps[-2] in current_dict:
                    current_dict[comps[-2]] += int(comps[-1])
                else:
                    current_dict[comps[-2]] = int(comps[-1])

    with open(os.path.join(out_dir, "stats_" + marker + ".txt"), 'w') as file:
        write_quality_information(quality_information, file)


def write_quality_information(quality_dict, file, prefix=[]):
    for key in sorted(quality_dict.keys(), key=lambda x: int(x)):
        if type(quality_dict[key]) == dict:
            write_quality_information(quality_dict[key], file, prefix + [key])
        else:
            file.write(SEPARATOR.join(prefix + [key, str(quality_dict[key])]) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Concatenates the statistics written by the alignment script into single files based on their suffixes (e.g. forward, reverse, total).')
    parser.add_argument('input', metavar='I', type=str,
                        help='The path to the statistics directory')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the output directory')
    args = parser.parse_args()

    markers = ['forward', 'reverse', 'total', 'length_deltas']
    for marker in markers:
        print("Concatenating {} files...".format(marker))
        combine_quality_stats(args.input, marker, args.output)
