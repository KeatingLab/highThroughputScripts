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

index_codes = ["ATCACG","CGATGT","TTAGGC","TGACCA","ACAGTG","GCCAAT",
                "CAGATC","ACTTGA","GATCAG"]
barcodes = ["ACTCG","ACTGT", "AATGC", "AGTCA", "ATACG", "ATAGC",
                "CGATC", "CTAAG", "CTCGA", "CGAAT", "CTGGT", "CGGTT",
                "GACTT", "GTTCA", "GATAC", "GAGCA", "GATGA", "GTCTG",
                "TCGGA", "TGACC", "TACTG", "TCCAG", "TCGAC", "TAGCT"]
FORMAT = "fastq-sanger"
SEQUENCE_QUALITY_KEY = "phred_quality"

BARCODE_FILE_PREFIX = "barcode_"

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

if __name__ == '__main__':
    # Split files. Decrease num_lines to increase the number of jobs into which the
    # task is split. The more jobs created, the more overhead needed to open file streams.
    print("Splitting files...")
    out_dir = "/Users/venkatesh-sivaraman/Documents/School/MIT/UROP2017/barcodes/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    splits_1 = os.path.join(out_dir, "split_1")
    splits_2 = os.path.join(out_dir, "split_2")
    test_1_paths = split_file("/Users/venkatesh-sivaraman/Documents/School/MIT/UROP2017/test_1.fastq", splits_1, num_lines=4000)
    test_2_paths = split_file("/Users/venkatesh-sivaraman/Documents/School/MIT/UROP2017/test_2.fastq", splits_2, num_lines=4000)

    #Process files
    print("Processing files...")
    path_pairs = zip(test_1_paths, test_2_paths)
    def processor(i):
        print("Processing chunk {}".format(i))
        test_1, test_2 = path_pairs[i]
        chunk_out_dir = os.path.join(out_dir, str(i))
        if not os.path.exists(chunk_out_dir):
            os.mkdir(chunk_out_dir)
        sort_sequence_reads(test_1, test_2, chunk_out_dir)

    pool = multiprocessing.Pool(processes=10)
    pool.imap_unordered(processor, xrange(len(path_pairs)))
    pool.close()
    pool.join()

    #Join files
    print("Joining files...")
    for i in xrange(len(index_codes) * len(barcodes)):
        #Get the paths for the appropriate reads for each chunk
        barcode_paths = [os.path.join(out_dir, str(chunk_number), BARCODE_FILE_PREFIX + str(i)) for chunk_number in xrange(len(path_pairs))]
        join_files(barcode_paths, os.path.join(out_dir, BARCODE_FILE_PREFIX + str(i)))

    if raw_input("Delete intermediate files? Type y for yes, n for no: ") == "y":
        paths_to_delete = [splits_1, splits_2] + [os.path.join(out_dir, str(i)) for i in xrange(len(path_pairs))]
        for path in paths_to_delete:
            if os.path.exists(path):
                print("Deleting {}".format(path))
                shutil.rmtree(path)
    print("Done.")
