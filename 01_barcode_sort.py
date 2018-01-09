'''
GOAL: Take .fastq files, corresponding to forward and reverse reads, and
sort each read into separate files depending on which barcodes it contains.
'''

from Bio import SeqIO
import sys
import os
import shutil
from itertools import izip
import multiprocessing
from functools import partial

index_codes = ["ATCACG","CGATGT","TTAGGC","TGACCA","ACAGTG","GCCAAT",
                "CAGATC","ACTTGA","GATCAG"]
barcodes = ["ACTCG","ACTGT", "AATGC", "AGTCA", "ATACG", "ATAGC",
                "CGATC", "CTAAG", "CTCGA", "CGAAT", "CTGGT", "CGGTT",
                "GACTT", "GTTCA", "GATAC", "GAGCA", "GATGA", "GTCTG",
                "TCGGA", "TGACC", "TACTG", "TCCAG", "TCGAC", "TAGCT"]
FORMAT = "fastq-sanger"
SEQUENCE_QUALITY_KEY = "phred_quality"

BARCODE_FILE_PREFIX = "barcode_"

NUM_PROCESSES = 15

def sort_sequence_reads(forward_path, reverse_path, out_dir):
    records_forward = SeqIO.parse(open(forward_path, "rU"), FORMAT)
    records_reverse = SeqIO.parse(open(reverse_path, "rU"), FORMAT)

    output_streams = open_output_streams(out_dir, BARCODE_FILE_PREFIX, len(index_codes) * len(barcodes))

    for (forward, reverse) in izip(records_forward, records_reverse):
        barcode_number = get_barcode_number(forward, reverse)
        if barcode_number == -1: continue

        SeqIO.write(forward, output_streams[barcode_number], FORMAT)
        SeqIO.write(reverse, output_streams[barcode_number], FORMAT)

    close_output_streams(output_streams)

    return 1

def open_output_streams(out_dir, prefix, num_files):
    return [open(os.path.join(out_dir, prefix + str(i)), "w") for i in xrange(num_files)]

def close_output_streams(streams):
    for stream in streams:
        stream.close()

def get_barcode_number(forward, reverse):
    '''
    Determines the file number for the given forward and reverse sequences.
    '''
    sequence = str(forward.seq)

    forward_quality = forward.letter_annotations[SEQUENCE_QUALITY_KEY]
    # If any element of the sequence has quality less than 20, discard the read
    if any(True for i in forward_quality if i < 20):
       return -1

    candidate_index = reverse.seq[:6].reverse_complement()
    # For now, demands an exact match
    index = get_closest_barcode(candidate_index, index_codes, 6)
    if index == -1:
       return -1

    candidate_barcode = sequence[:5]
    # Exact match
    barcode = get_closest_barcode(candidate_barcode, barcodes, 5)
    if barcode == -1:
       return -1

    return (index * len(barcodes)) + barcode

def accuracy(seq, other_seq, threshold):
    '''
    Returns the number of identical nucleotides between seq and other_seq. If the
    accuracy is less than threshold, returns -1.
    '''
    acc = len([i for i in range(len(seq)) if seq[i] == other_seq[i]])
    if acc < threshold:
        return -1
    return acc

def get_closest_barcode(sequence, barcode_list, match_threshold=None):
    '''
    Returns the index of the element of barcode_list that most closely matches
    sequence. The number of identical nucleotides must be at least match_threshold.
    If match_threshold is not passed, then any number of identical nucleotides
    will be admitted.
    '''
    threshold = match_threshold if match_threshold is not None else 0
    return max(xrange(len(barcode_list)), key=lambda i: accuracy(sequence, barcode_list[i], threshold))

#### Splitting files

def split_file(path, out_dir, num_lines=2e7):
    '''
    Splits the given file into smaller files in the given output directory, where
    each smaller file has at most the given number of lines. Returns the list of
    file paths created.
    '''
    file_name, extension = os.path.splitext(path)
    file_name = os.path.basename(file_name)
    ret = []

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    with open(path, "rU") as file:
        current_out_file = None
        current_out_index = 0
        line_count = 0

        for line in file:
            if current_out_file is None:
                new_path = os.path.join(out_dir, file_name + "_" + str(current_out_index) + extension)
                current_out_file = open(new_path, "w")
                current_out_index += 1
                ret.append(new_path)

            current_out_file.write(line)
            line_count += 1

            if line_count >= num_lines:
                current_out_file.close()
                current_out_file = None
                line_count = 0

        if current_out_file is not None:
            current_out_file.close()

    return ret

def join_files(file_paths, out_file):
    '''
    Concatenates the files at the given paths into the file at out_file.
    '''
    with open(out_file, "w") as file:
        for path in file_paths:
            with open(path, "rU") as chunk:
                for line in chunk:
                    file.write(line)

def sort_barcodes_processor(out_dir, input_info):
    '''
    Worker function for the multiprocessing pool.
    '''
    i, test_1, test_2 = input_info
    print("Processing chunk {}".format(i))
    chunk_out_dir = os.path.join(out_dir, str(i))
    if not os.path.exists(chunk_out_dir):
        os.mkdir(chunk_out_dir)
    sort_sequence_reads(test_1, test_2, chunk_out_dir)

def split_files_processor(num_lines, input_info):
    path, out_dir = input_info
    return split_file(path, out_dir, num_lines)

def join_files_processor(out_dir, num_chunks, barcode):
    #Get the paths for the appropriate reads for each chunk
    barcode_paths = [os.path.join(out_dir, str(chunk_number), BARCODE_FILE_PREFIX + str(barcode)) for chunk_number in xrange(num_chunks)]
    join_files(barcode_paths, os.path.join(out_dir, BARCODE_FILE_PREFIX + str(barcode)))


def sort_barcodes_with_split(forward_path, reverse_path, out_dir):
    # Split files. Decrease num_lines to increase the number of jobs into which the
    # task is split. The more jobs created, the more overhead needed to open file streams.
    print("Splitting files...")
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    splits_1 = os.path.join(out_dir, "split_1")
    splits_2 = os.path.join(out_dir, "split_2")
    pool = multiprocessing.Pool(processes=2)
    processor = partial(split_files_processor, 4000) # Number of lines per file - must be 0 mod 4
    #Input, output pairs for each splitting job
    split_args = [(forward_path, splits_1), (reverse_path, splits_2)]
    forward_paths, reverse_paths = pool.map(processor, split_args)

    #Process files
    input_items = [(i, forward_paths[i], reverse_paths[i]) for i in range(len(forward_paths))]

    print("Processing files ({} chunks)...".format(len(input_items)))
    pool = multiprocessing.Pool(processes=NUM_PROCESSES)
    processor = partial(sort_barcodes_processor, out_dir)
    # Synchronous equivalent:
    # map(processor, input_items)
    result = pool.imap(processor, input_items)
    pool.close()
    pool.join()

    print("Joining files...")
    pool = multiprocessing.Pool(processes=NUM_PROCESSES)
    processor = partial(join_files_processor, out_dir, len(input_items))
    pool.map(processor, range(len(index_codes) * len(barcodes)))

    print("Cleaning up...")
    paths_to_delete = [splits_1, splits_2] + [os.path.join(out_dir, str(i)) for i in xrange(len(input_items))]
    for path in paths_to_delete:
        if os.path.exists(path):
            print("Deleting {}".format(path))
            shutil.rmtree(path)
    print("Done.")


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Not enough arguments. Provide the path to the forward reads, then the path to the reverse reads, then the path to the output directory.")
        exit(1)

    in_path_1 = sys.argv[1]
    in_path_2 = sys.argv[2]
    out_dir = sys.argv[3]
    sort_barcodes_with_split(in_path_1, in_path_2, out_dir)
