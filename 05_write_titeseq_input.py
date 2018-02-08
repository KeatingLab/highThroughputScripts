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

TASKS = [
    (TASK_TRACK_SINGLE_SEQ, ["barcode_0_0", "barcode_1_0"], "AAACAAGAACCTCAGGAAATCGATTTCCCGGACGATCTGCC", "wt_counts.csv"),
    (TASK_TRACK_SEQS, ["barcode_0_0", "barcode_1_0"], "titeseq_input_1.csv")
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
        if return_counts:
            total_reads += count
        if target_sequence is not None and root != target_sequence:
            continue
        if root not in processed_counts:
            processed_counts[root] = {}
        processed_counts[root][basename] = count

    if return_counts:
        return processed_counts, len(seqs), total_reads
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
        for task in TASKS:
            perform_task(task, args)
    else:
        task = TASKS[args.task]
        perform_task(task, args)

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))

    # Write out the input parameters to file
    params_path = os.path.join(args.output, "params.txt")
    sc.write_input_parameters([(name, eval(name)) for name in PARAMETER_LIST], params_path)
