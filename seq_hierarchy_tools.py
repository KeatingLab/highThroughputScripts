INDENT_CHARACTER = '\t'
DELIMITER_CHARACTER = '\t'
INDEX_HELPER_KEY = "__index__"

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
