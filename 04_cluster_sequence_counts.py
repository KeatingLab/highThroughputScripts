'''
This script operates on the complete output of 03_count_sequences, which is
a list of frequencies of unique sequences subindexed by the frequencies of other
unique sequences within the same reads. It will first cluster the unique top-level
sequences by a similarity cutoff.
'''

import os, sys
import time
import argparse
from aligner import Aligner
import multiprocessing
from functools import partial

INDENT_CHARACTER = '\t'
DELIMITER_CHARACTER = '\t'

'''
Each tuple should contain the following:
* the mode of operation ("cluster")
* the hierarchical level at which to perform the operation (0 is top-level)
* the number of mutations that should be allowed
* a list of ranges in which the sequences should be compared, or None
'''
TASKS = [("cluster", 0, 1, [(0, 27)]),
         ("cluster", 1, 1, None)]

def read_sequence_dicts(in_file):
    '''
    Reads all of the sequences and counts from the given file and returns them as
    a list of dictionaries. For instance, if 03_count_sequences.py were used to
    group a set of sequences by Region A, then by Region B, the output of this
    method would look like the following:
    [
    { A1: [ { B1: 3 }, { B2: 4 }, { B3: 1 }, ...] },
    { A2: [ { B4: 5 }, { B5: 3 }, { B6: 2 }, ...] },
    ...
    ]
    '''
    ret = []
    current_reading_sequence = {}
    current_key_path = []
    current_indent_level = 0
    for line in in_file:
        indent_level = len(line) - len(line.lstrip(INDENT_CHARACTER))
        count, unique_count, seq = line.strip().split(DELIMITER_CHARACTER)
        if indent_level == 0:
            # Add the currently-built entry and initialize a new empty one
            if len(current_reading_sequence) > 0:
                ret.append(current_reading_sequence)
            current_reading_sequence = {seq: int(count)}
            current_key_path = [seq]
        else:
            if indent_level > current_indent_level:
                # Create a new sub-dictionary
                set_at_key_path(current_reading_sequence, current_key_path, [{seq: int(count)}])
                current_key_path += [0, seq]
            else:
                if indent_level < current_indent_level:
                    # Move into the outer dictionary
                    current_key_path.pop()
                    current_key_path.pop()
                else:
                    # Create a sibling element
                    current_key_path[-2] += 1
                    current_key_path[-1] = seq
                current_item = get_at_key_path(current_reading_sequence, current_key_path[:-2])
                current_item.append({seq: int(count)})
        current_indent_level = indent_level

    if len(current_reading_sequence) > 0:
        ret.append(current_reading_sequence)
    return ret

def write_sequence_dicts(sequences, out_path, out_file=None, indent_level=0):
    '''
    Writes the given list of sequence hierarchy dictionaries to the given file.
    Uses the same format as the output of 03_count_sequences.py.
    '''
    if out_file is None:
        file = open(out_path, 'w')
    else:
        file = out_file

    for sequence_dict in sequences:
        seq = sequence_dict.keys()[0]
        value = sequence_dict[seq]
        try:
            count = int(value)
            file.write(INDENT_CHARACTER * (indent_level * 2) + DELIMITER_CHARACTER.join([str(count), str(count), seq]) + '\n')
        except:
            count, unique = get_sequence_counts(value)
            file.write(INDENT_CHARACTER * (indent_level * 2) + DELIMITER_CHARACTER.join([str(count), str(unique), seq]) + '\n')
            write_sequence_dicts(value, out_path, out_file, indent_level + 1)

    if out_file is None:
        file.close()

def get_sequence_counts(sequence_info):
    '''
    Returns the total count and the unique count given a list from a sequence
    info object.
    '''
    unique = len(sequence_info)
    total = 0
    for item in sequence_info:
        key, value = item.items()[0]
        try:
            total += int(value)
        except:
            total += get_sequence_counts(value)[0]
    return total, unique

def get_at_key_path(dictionary, key_path):
    '''
    Returns the element of the given nested dictionary at the given key path (a
    list of keys used to drill down into the dictionary).
    '''
    current = dictionary
    for key in key_path:
        current = current[key]
    return current

def set_at_key_path(dictionary, key_path, value):
    '''
    Sets the element of the given nested dictionary at the given key path.
    '''
    current = dictionary
    for key in key_path[:-1]:
        current = current[key]
    current[key_path[-1]] = value

def merge_sequence_info(source, new):
    '''
    Merges the two sequence items (counts or lists of dictionaries) and returns the result.
    '''
    try:
        return int(source) + int(new)
    except:
        pass
    ret = {}
    for item in source:
        key, value = item.items()[0]
        ret[key] = value
    for item in new:
        key, value = item.items()[0]
        if key in ret:
            if isinstance(ret[key], list):
                assert isinstance(value, list), "Mismatched types while merging: {} and {}".format(ret[key], value)
                ret[key] = merge_sequence_info(ret[key], value)
            else:
                ret[key] += value
        else:
            ret[key] = value
    return [{key: value} for key, value in ret.items()]

def compare_sequence_dictionaries(source, new, task):
    '''
    Returns True if the two sequence dictionaries should be merged.
    '''
    aligner = Aligner(different_score=0)
    source_key = source.keys()[0]
    new_key = new.keys()[0]
    scoring_maps = None
    max_score = len(source_key)
    _, _, allowed_mutations, scoring_ranges = task
    if scoring_ranges is not None:
        score_map = aligner.scoring_map(len(source_key), scoring_ranges)
        scoring_maps = (score_map, score_map)
        max_score = sum(score_map)
    return aligner.score(source_key, new_key, 0, scoring_maps=scoring_maps) >= max_score - allowed_mutations

def cluster_sequences_processor(task, seq, (index, other_seq)):
    '''
    Convenience function for multiprocessing that calls the merge_function and
    returns part of the input.
    '''
    return index, compare_sequence_dictionaries(seq, other_seq, task)

def cluster_sequences(all_sequences, task, num_processes=15, current_level=0):
    '''
    Returns an iterator that, for each pair of candidate sequences, calls the
    merge_function with two arguments (source_seqs,
    candidate_seqs), where each is a nested dictionary keyed by sequences
    referenced in the input file and where the values are either sub-
    dictionaries or counts. If the merge_function returns None, the candidate
    sequences are kept separate from the source sequences. If the function
    returns another sequence dictionary, this dictionary will be used as the
    merged result and yielded.
    '''
    visited_indexes = set()
    pool = multiprocessing.Pool(processes=num_processes)

    level = task[1]

    for i in xrange(len(all_sequences)):
        if i in visited_indexes:
            continue

        seq = all_sequences[i]
        root = seq.keys()[0]
        visited_indexes.add(i)

        if level == current_level:
            merged_seq = seq
            processor = partial(cluster_sequences_processor, task, seq)
            indexes = ((j, all_sequences[j]) for j in xrange(i + 1, len((all_sequences))) if j not in visited_indexes)
            for index, result in map(processor, indexes): #pool.imap(processor, indexes, chunksize=1000):
                if result:
                    other_seq = all_sequences[index]
                    other_root = other_seq.keys()[0]
                    merged_seq[root] = merge_sequence_info(merged_seq[root], other_seq[other_root])
                    visited_indexes.add(index)
            print("After {}, visited {} sequences".format(i, len(visited_indexes)))

        elif level > current_level:
            merged_seq = {root: []}
            for result in cluster_sequences(seq[root], task, num_processes, current_level + 1):
                merged_seq[root].append(result)

        yield merged_seq
        if len(visited_indexes) == len(all_sequences):
            break

def cluster_sequences_from_file(in_path, out_dir, tasks, **kwargs):
    '''
    Performs cluster_sequences on the contents of the file at in_path. The
    resulting merged data is written to the appropriate files within out_dir.

    This is a worst-case quadratic-time algorithm and loads all sequences in the
    input file into memory.
    '''
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for i, task in enumerate(tasks):
        print("Task {}".format(i))
        with open(in_path, 'r') as in_file:
            all_sequences = read_sequence_dicts(in_file)

        task_dir = os.path.join(out_dir, "task_{}".format(i))
        if not os.path.exists(task_dir):
            os.mkdir(task_dir)
        out_path = os.path.join(task_dir, os.path.basename(in_path))

        with open(out_path, 'w') as out_file:
            for cluster in cluster_sequences(all_sequences, task, **kwargs):
                write_sequence_dicts([cluster], out_path, out_file)
        # Use out_path as input for next stage
        in_path = out_path


if __name__ == '__main__':
    a = time.time()  # Time the script started
    parser = argparse.ArgumentParser(description='Clusters the given sequence count hierarchy (e.g. the output of 03_count_sequences.py) according to the tasks listed at the top of the script. One file will be written out for each task.')
    parser.add_argument('input', metavar='I', type=str,
                        help='The path to the input sequence hierarchy file')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the output directory')
    parser.add_argument('-p', '--processes', type=int, default=15,
                        help='The number of processes to use')
    args = parser.parse_args()

    cluster_sequences_from_file(args.input, args.output, TASKS, num_processes=args.processes)

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))
