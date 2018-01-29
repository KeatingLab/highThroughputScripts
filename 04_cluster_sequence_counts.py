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
INDEX_HELPER_KEY = "__index__"

TASK_CLUSTER = "cluster"
TASK_SWITCH = "switch"
TASK_PASS = "pass"

'''
Each tuple should contain the following:
* the mode of operation
* if the mode is TASK_CLUSTER:
    * the hierarchical level at which to perform the operation (0 is top-level)
    * the number of mutations that should be allowed
    * a list of ranges in which the sequences should be compared, or None
* if the mode is TASK_SWITCH:
    * the hierarchical level which should be moved to the top level
* if the mode is TASK_PASS:
    * no additional arguments. Make sure the item is still a tuple - e.g., use
        `(TASK_PASS,)`. This is a useful parameter to pass if you want to proceed
        forward using previously-generated intermediate output.
'''
TASKS = [(TASK_CLUSTER, 0, 1, [(0, 27)]),
         (TASK_CLUSTER, 1, 1, None),
         (TASK_SWITCH, 1)]
         #(TASK_CLUSTER, 0, 1, None),
         #(TASK_SWITCH, 1)]

'''
List of expression strings which will be evaluated and written out to the
params.txt file.
'''
PARAMETER_LIST = ["args.input", "args.output", "args.processes", "TASKS"]

### I/O for sequence dictionary objects

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
        seq = get_root_item(sequence_dict)[0]
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

def get_at_key_path(dictionary, key_path):
    '''
    Returns the element of the given nested dictionary at the given key path (a
    list of keys used to drill down into the dictionary). Used for reading and
    writing sequence dictionaries.
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

### Helper functions

def get_root_item(sequence):
    '''
    Returns the root key and value in the given dictionary.
    '''
    for key, value in sequence.iteritems():
        if key != INDEX_HELPER_KEY:
            return key, value
    return None

def get_sequence_counts(sequence_info):
    '''
    Returns the total count and the unique count given a list from a sequence
    info object.
    '''
    unique = len(sequence_info)
    total = 0
    for item in sequence_info:
        key, value = get_root_item(item)
        try:
            total += int(value)
        except:
            total += get_sequence_counts(value)[0]
    return total, unique

def get_children_at_level(sequence, level, current_level=0):
    '''
    Yields each sequence hierarchy object (dictionary) at the given level along
    with the key path that led to it. Used for performing switches of the
    hierarchy level.
    '''
    if current_level == level:
        yield [], sequence
    else:
        root_key, root_value = get_root_item(sequence)
        if type(root_value) is not int:
            for subseq in root_value:
                for key_path, child in get_children_at_level(subseq, level, current_level + 1):
                    yield [root_key] + key_path, child

def set_tree_at_key_path(sequence, key_path, tree):
    '''
    Given a list of keys in key_path, traverse down into the sequence dictionary
    and add the given tree at that key path.
    '''
    current_dict = sequence
    assert len(key_path) >= 1, "Cannot set a count without at least one key"
    for i, key in enumerate(key_path):
        root_key, root_value = get_root_item(current_dict)
        if INDEX_HELPER_KEY not in current_dict:
            current_dict[INDEX_HELPER_KEY] = {}
            for i, item in enumerate(root_value):
                sub_root_key = get_root_item(item)[0]
                current_dict[INDEX_HELPER_KEY][sub_root_key] = i
            #print("Initialized index for {}".format(current_dict))

        if key in current_dict[INDEX_HELPER_KEY]:
            current_dict = root_value[current_dict[INDEX_HELPER_KEY][key]]
        else:
            current_dict[INDEX_HELPER_KEY][key] = len(root_value)
            current_dict = {key: []}
            root_value.append(current_dict)

        if i == len(key_path) - 1:
            current_dict[key] = tree

### Clustering and merging

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
        key, value = get_root_item(item)
        ret[key] = value
    for item in new:
        key, value = get_root_item(item)
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
    source_key = get_root_item(source)[0]
    new_key = get_root_item(new)[0]
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
    level = task[1]

    if level == current_level and len(all_sequences) == 1:
        yield all_sequences[0]
        return

    visited_indexes = set()
    pool = multiprocessing.Pool(processes=num_processes)

    for i in xrange(len(all_sequences)):
        if i in visited_indexes:
            continue

        seq = all_sequences[i]
        root = get_root_item(seq)[0]
        visited_indexes.add(i)

        if level == current_level:
            merged_seq = seq
            processor = partial(cluster_sequences_processor, task, seq)
            indexes = ((j, all_sequences[j]) for j in xrange(i + 1, len((all_sequences))) if j not in visited_indexes)
            for index, result in map(processor, indexes):#pool.imap(processor, indexes, chunksize=1000):
                if result:
                    other_seq = all_sequences[index]
                    other_root = get_root_item(other_seq)[0]
                    merged_seq[root] = merge_sequence_info(merged_seq[root], other_seq[other_root])
                    visited_indexes.add(index)
            if level == 0:
                print("After {}, visited {} sequences".format(i, len(visited_indexes)))

        elif level > current_level:
            merged_seq = {root: []}
            for result in cluster_sequences(seq[root], task, num_processes, current_level + 1):
                merged_seq[root].append(result)

        yield merged_seq
        if len(visited_indexes) == len(all_sequences):
            break

    pool.close()
    pool.join()

### Switching sequence hierarchy

def switch_sequence_hierarchy(all_sequences, task):
    '''
    Generates a new sequence hierarchy (nested dictionaries and lists), where the
    level specified by the given task has been moved to the top.
    '''
    new_hierarchy = []
    top_level_indexes = {}
    level = task[1]
    for sequence in all_sequences:
        for key_path, child in get_children_at_level(sequence, level):
            # Example:
            # key_path: ['TGGTGTCATCCTTGTGAAAGTGATTTTCCGGACGATCTGCC', 'GGTATATAA']
            # child: {'ACTCG': 1}
            root_key, root_value = get_root_item(child)
            if root_key not in top_level_indexes:
                top_level_indexes[root_key] = len(new_hierarchy)
                new_hierarchy.append({root_key: []})

            index_of_key = top_level_indexes[root_key]
            set_tree_at_key_path(new_hierarchy[index_of_key], key_path, root_value)
    for cluster in new_hierarchy:
        #print(cluster)
        yield cluster

### Main

def collapse_sequences_from_file(in_path, out_dir, tasks, **kwargs):
    '''
    Performs the given set of tasks on the sequence hierarchies contained at
    in_path. Each successive task is performed on the output of the previous task.

    This function uses worst-case quadratic-time algorithms and loads all sequences in the
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

        if task[0] != TASK_PASS:
            with open(out_path, 'w') as out_file:
                if task[0] == TASK_CLUSTER:
                    for cluster in cluster_sequences(all_sequences, task, **kwargs):
                        write_sequence_dicts([cluster], out_path, out_file)
                elif task[0] == TASK_SWITCH:
                    for cluster in switch_sequence_hierarchy(all_sequences, task):
                        write_sequence_dicts([cluster], out_path, out_file)
                else:
                    print("Undefined task: {}".format(task))
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

    collapse_sequences_from_file(args.input, args.output, TASKS, num_processes=args.processes)

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))

    # Write out the input parameters to file
    params_path = os.path.join(args.output, "params.txt")
    sc.write_input_parameters([(name, eval(name)) for name in PARAMETER_LIST], params_path)
