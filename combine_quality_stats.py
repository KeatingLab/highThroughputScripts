'''
Concatenates CSV files such as those created by the stat_collector script.
'''

import os
import argparse

SEPARATOR = '\t'

def combine_quality_stats(input_dir, marker, out_dir):
    quality_information = {}
    file_count = 0
    for path in os.listdir(input_dir):
        if path[0] == '.': continue
        if len(marker) > 0 and marker not in path: continue
        file_count += 1
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
    print("Concatenated {} files with marker '{}'.".format(file_count, marker))


def write_quality_information(quality_dict, file, prefix=[]):
    for key in sorted(quality_dict.keys(), key=lambda x: int(x)):
        if type(quality_dict[key]) == dict:
            write_quality_information(quality_dict[key], file, prefix + [key])
        else:
            file.write(SEPARATOR.join(prefix + [key, str(quality_dict[key])]) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Concatenates statistics written by the stat_collector module into single files. Can either concatenate the entire input directory, or split by markers in the input file names (e.g. 'forward', 'reverse', etc.).")
    parser.add_argument('input', metavar='I', type=str,
                        help='The path to the statistics directory')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the output directory')
    parser.add_argument('markers', metavar='M', type=str, nargs='*', default=[""],
                        help='The markers for each type of statistic to use (if nothing is passed, concatenates all contents of the input directory)')
    args = parser.parse_args()

    for marker in args.markers:
        combine_quality_stats(args.input, marker, args.output)
