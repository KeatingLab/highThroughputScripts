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

TASKS = [
    (TASK_TRACK_SEQS, ["barcode_0_0", "barcode_1_0"], "titeseq_input_1.csv")
]

'''
List of expression strings which will be evaluated and written out to the
params.txt file.
'''
PARAMETER_LIST = ["args.input", "args.output", "TASKS"]

def process_sequence_file(in_path, processed_counts):
    '''
    in_path is a path to a sequence counts file, and processed_counts is a
    dictionary of sequences mapped to dictionaries of file names to counts. For
    instance:
    {'ABCDE': {'sequence_file_1': 3, 'sequence_file_2': 5, ...},
     'BCDEF': {'sequence_file_1': 1, 'sequence_file_2': 7, ...}}
    Returns the updated processed_counts dictionary.
    '''
    basename = os.path.basename(in_path)

    with open(in_path, 'r') as file:
        seqs = read_sequence_dicts(file)

    for seq in seqs:
        root, value = get_root_item(seq)
        count, unique = get_sequence_counts(value)
        if root not in processed_counts:
            processed_counts[root] = {}
        processed_counts[root][basename] = unique

    return processed_counts

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

    with open(out_path, "w") as file:
        for seq, counts in result.items():
            comps_list = [seq]
            for path in paths:
                if path in counts:
                    comps_list.append(str(counts[path]))
                else:
                    comps_list.append("0")
            file.write(",".join(comps_list) + "\n")

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
        for task in TASKS:
            if task[0] == TASK_TRACK_SEQS:
                track_sequences(args.input, args.output, task)
            else:
                print("Unrecognized task type {}".format(task[0]))
    else:
        task = TASKS[args.task]
        if task[0] == TASK_TRACK_SEQS:
            track_sequences(args.input, args.output, task)
        else:
            print("Unrecognized task type {}".format(task[0]))

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))

    # Write out the input parameters to file
    params_path = os.path.join(args.output, "params.txt")
    sc.write_input_parameters([(name, eval(name)) for name in PARAMETER_LIST], params_path)
