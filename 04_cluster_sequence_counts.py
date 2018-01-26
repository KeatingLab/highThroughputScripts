'''
This script operates on the complete output of 03_count_sequences, which is
a list of frequencies of unique sequences subindexed by the frequencies of other
unique sequences within the same reads. It will first cluster the unique top-level
sequences by a similarity cutoff.
'''

import sys
from aligner import Aligner

INDENT_CHARACTER = '\t'
DELIMITER_CHARACTER = '\t'

ALLOWED_SEQUENCE_MUTATIONS = 1

def read_sequence_dicts(in_file):
    '''
    Reads all of the sequences and counts from the given file and returns them as
    a list of dictionaries. For instance, if 03_count_sequences.py were used to
    group a set of sequences by Region A, then by Region B, the output of this
    method would look like the following:
    [
    { A1: { B1: 3, B2: 4, B3: 1, ...} },
    { A2: { B4: 5, B5: 3, B6: 2, ...} },
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
                set_at_key_path(current_reading_sequence, current_key_path, {seq: int(count)})
                current_key_path.append(seq)
            else:
                if indent_level < current_indent_level:
                    # Move into the outer dictionary
                    current_key_path.pop()
                else:
                    # Create a sibling element
                    current_key_path[-1] = seq
                current_item = get_at_key_path(current_reading_sequence, current_key_path[:-1])
                current_item[seq] = int(count)
                set_at_key_path(current_reading_sequence, current_key_path[:-1], current_item)
        current_indent_level = indent_level

    if len(current_reading_sequence) > 0:
        ret.append(current_reading_sequence)
    return ret

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

def cluster_sequences(all_sequences, merge_function):
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
    for i in xrange(len(all_sequences)):
        if i in visited_indexes:
            continue
        seq = all_sequences[i]
        root = seq.keys()[0]
        merged_seq = seq
        visited_indexes.add(i)

        for j in xrange(i + 1, len((all_sequences))):
            if j in visited_indexes:
                continue
            other_seq = all_sequences[j]
            other_root = other_seq.keys()[0]
            result = merge_function(merged_seq, other_seq)
            if result is not None:
                merged_seq = result
                visited_indexes.add(j)
        yield merged_seq
        print("After {}, visited {} sequences".format(i, len(visited_indexes)))
        if len(visited_indexes) == len(all_sequences):
            break

def cluster_sequences_from_file(in_file, out_file, merge_function):
    '''
    Performs cluster_sequences on the contents of in_file. The resulting merged
    data is written to the out_file.

    This is a worst-case quadratic-time algorithm and loads all sequences in the
    input file into memory.
    '''
    all_sequences = read_sequence_dicts(in_file)
    merged_sequences = []
    for cluster in cluster_sequences(all_sequences, merge_function):
        merged_sequences.append(cluster)
    print("Original length: {} now: {}".format(len(all_sequences), len(merged_sequences)))
    print('\n'.join(str(s) for s in merged_sequences))

def merge_sequence_dicts(source, new):
    '''
    Merges the two sequence dictionaries and returns the result.
    '''
    ret = {}
    for key, value in source.items():
        if key in new:
            if type(value) is dict:
                assert type(new[key]) is dict, "Mismatched types between sequence dictionaries: {} and {}".format(type(value), type(new[key]))
                ret[key] = merge_sequence_dicts(value, new[key])
            else:
                ret[key] = value + new[key]
        else:
            ret[key] = value
    for key, value in new.items():
        if key not in source:
            ret[key] = value
    return ret

def compare_sequence_dictionaries(source, new):
    '''
    Returns a merged dictionary if the two sequence dictionaries should be merged.
    '''
    aligner = Aligner(different_score=0)
    source_key = source.keys()[0]
    new_key = new.keys()[0]
    if aligner.score(source_key, new_key, 0) >= len(source_key) - ALLOWED_SEQUENCE_MUTATIONS:
        return {source_key: merge_sequence_dicts(source[source_key], new[new_key])}
    return None

if __name__ == '__main__':
    input_path = sys.argv[1]
    with open(input_path, 'r') as file:
        cluster_sequences_from_file(file, None, compare_sequence_dictionaries)
