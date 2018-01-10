'''
GOAL: Take .fastq files, corresponding to forward and reverse reads, and
sort each read into separate files depending on which barcodes it contains.
'''

from Bio import SeqIO
import sys
import os
import shutil
from itertools import izip
import multiprocessing, threading
from functools import partial
import time

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

def sort_sequence_reads(forward_path, reverse_path, out_dir, append=False, threads=1, chunk_size=1):
    records_forward = SeqIO.parse(open(forward_path, "rU"), FORMAT)
    records_reverse = SeqIO.parse(open(reverse_path, "rU"), FORMAT)

    output_streams = open_output_streams(out_dir, BARCODE_FILE_PREFIX, len(index_codes) * len(barcodes), append)

    if threads > 1:
        manager = multiprocessing.Manager()
        processor = partial(sort_sequence, output_streams)
        pool = multiprocessing.Pool(processes=threads)
        pool.imap(processor, enumerate(izip(records_forward, records_reverse)))
        pool.close()
        pool.join()
    elif chunk_size > 1:
        sort_sequence_chunks(records_forward, records_reverse, chunk_size, output_streams)
    else:
        processor = partial(sort_sequence, output_streams)
        map(processor, enumerate(izip(records_forward, records_reverse)))

    close_output_streams(output_streams)

    return 1

def sort_sequence(output_streams, (index, (forward, reverse)), locks=None):
    if index % 100 == 0:
        print("Processing item {}".format(index))
    barcode_number = get_barcode_number(forward, reverse)
    if barcode_number == -1: return

    if locks is not None:
        locks[barcode_number].acquire()
    SeqIO.write(forward, output_streams[barcode_number], FORMAT)
    SeqIO.write(reverse, output_streams[barcode_number], FORMAT)
    if locks is not None:
        locks[barcode_number].release()

def sort_sequence_chunks(records_forward, records_reverse, chunk_size, output_streams):
    collected_forward = []
    def collate_records():
        # Collect reverse records now
        barcode_bins = {}
        for record_index in xrange(len(collected_forward)):
            reverse = next(records_reverse, None)
            if reverse is None: break
            barcode_number = get_barcode_number(collected_forward[record_index], reverse)
            if barcode_number == -1: continue
            if barcode_number in barcode_bins:
                barcode_bins[barcode_number].append((forward, reverse))
            else:
                barcode_bins[barcode_number] = [(forward, reverse)]
        # Write the newly completed entries to disk
        for barcode, items in barcode_bins.items():
            for forward_item, reverse_item in items:
                SeqIO.write(forward_item, output_streams[barcode], FORMAT)
                SeqIO.write(reverse_item, output_streams[barcode], FORMAT)

    chunk_index = 0
    for forward in records_forward:
        collected_forward.append(forward)
        if len(collected_forward) == chunk_size:
            print("Writing chunk {} to file...".format(chunk_index))
            collate_records()
            collected_forward = []
            chunk_index += 1
    collate_records()

def open_output_streams(out_dir, prefix, num_files, append=False):
    return [open(os.path.join(out_dir, prefix + str(i)), "a" if append else "w") for i in xrange(num_files)]

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

def write_file_chunk(file, output, num_lines):
    '''
    Returns whether or not the file stream ended with this chunk.
    '''
    dest = os.path.dirname(output)
    if not os.path.exists(dest):
        os.mkdir(dest)
    line_count = 0
    with open(output, "w") as out_file:
        while line_count < num_lines:
            line = file.readline()
            if len(line) == 0:
                return True
            out_file.write(line)
            line_count += 1
    return False

def join_files(file_paths, out_file):
    '''
    Concatenates the files at the given paths into the file at out_file.
    '''
    with open(out_file, "w") as file:
        for path in file_paths:
            if not os.path.exists(path): continue
            with open(path, "rU") as chunk:
                for line in chunk:
                    file.write(line)

def join_files_processor(out_dir, num_chunks, barcode):
    #Get the paths for the appropriate reads for each chunk
    barcode_paths = [os.path.join(out_dir, str(chunk_number), BARCODE_FILE_PREFIX + str(barcode)) for chunk_number in xrange(num_chunks)]
    join_files(barcode_paths, os.path.join(out_dir, BARCODE_FILE_PREFIX + str(barcode)))
    for path in barcode_paths:
        if not os.path.exists(path): continue
        os.remove(path)

class StreamingSortWorker(object):
    '''
    An object that lives in its own thread and sorts the reads from the fastq
    files it is given into barcodes in its own out directory.
    '''
    def __init__(self, id, out_base_dir):
        self.id = id
        self.out_dir = os.path.join(out_base_dir, str(id))
        self.is_working = False
        self.chunk_number = -1
        self.thread = None

    def sort(self, forward_input, reverse_input, delete_on_complete=True):
        self.thread = threading.current_thread()
        self.is_working = True
        if not os.path.exists(self.out_dir):
            os.mkdir(self.out_dir)
        sort_sequence_reads(forward_input, reverse_input, self.out_dir, append=True)
        if delete_on_complete:
            os.remove(forward_input)
            os.remove(reverse_input)
        self.is_working = False
        self.thread = None

available_chunk_files = []

def streaming_work_supplier(lock, streams, split_dirs, num_lines, num_chunk_files):
    global available_chunk_files
    # Write the chunks
    forward_stream, reverse_stream = streams
    chunk_number = 0
    while True:
        with lock:
            available_files = len(available_chunk_files)

        if available_files < num_chunk_files:
            print("Transferring chunk {}...".format(chunk_number))
            forward_path = os.path.join(split_dirs[0], str(chunk_number))
            reverse_path = os.path.join(split_dirs[1], str(chunk_number))
            finished_1 = write_file_chunk(forward_stream, forward_path, num_lines)
            finished_2 = write_file_chunk(reverse_stream, reverse_path, num_lines)
            with lock:
                available_chunk_files.append((forward_path, reverse_path))

            chunk_number += 1
            if finished_1 or finished_2:
                print("Finished reading chunks.")
                return

def sort_barcodes_streaming(forward_path, reverse_path, out_dir, num_lines=6e7, num_threads=4):
    '''
    The idea behind the streaming version is to transfer chunks of size num_lines
    into temporary files, then hand them to worker threads to sort. Each worker
    thread handles its own set of output files, which will be joined together at
    the end of this function. The "master" (this method) is responsible for making
    sure that every worker always has a chunk to work on, and that no chunk is
    being handled by multiple workers.
    '''
    global available_chunk_files

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    splits = (os.path.join(out_dir, "split_1"), os.path.join(out_dir, "split_2"))
    lock = threading.Lock()
    streams = (open(forward_path, "rU"), open(reverse_path, "rU"))

    workers = [StreamingSortWorker(i, out_dir) for i in xrange(num_threads)]
    supply_thread = threading.Thread(target=streaming_work_supplier, args=(lock, streams, splits, num_lines, num_threads * 3))
    supply_thread.start()

    while supply_thread.is_alive() or len(available_chunk_files) > 0:
        for worker in workers:
            if not worker.is_working:
                with lock:
                    if len(available_chunk_files) == 0: continue
                    forward, reverse = available_chunk_files.pop(0)

                # Spawn a worker thread
                print("Starting chunk {}...".format(os.path.basename(forward)))
                worker_thread = threading.Thread(target=worker.sort, args=(forward, reverse))
                worker_thread.start()
        time.sleep(1)

    for worker in workers:
        if worker.thread is not None:
            worker.thread.join()
    map(file.close, streams)

    print("Joining files...")
    pool = multiprocessing.Pool(processes=NUM_PROCESSES)
    processor = partial(join_files_processor, out_dir, len(workers))
    pool.map(processor, range(len(index_codes) * len(barcodes)))

    print("Cleaning up...")
    paths_to_delete = list(splits) + [os.path.join(out_dir, str(i)) for i in xrange(len(workers))]
    for path in paths_to_delete:
        if os.path.exists(path):
            print("Deleting {}".format(path))
            shutil.rmtree(path)
    print("Done.")


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Not enough arguments. Provide the path to the forward reads, then the path to the reverse reads, then the path to the output directory.")
        exit(1)

    a = time.time()  # Time the script started
    in_path_1 = sys.argv[1]
    in_path_2 = sys.argv[2]
    out_dir = sys.argv[3]
    if len(sys.argv) > 4:
        num_lines = int(sys.argv[4])
    else:
        num_lines = 2e7
    #sort_barcodes_streaming(in_path_1, in_path_2, out_dir, num_lines=40000)
    sort_sequence_reads(in_path_1, in_path_2, out_dir, chunk_size=20000)
    b = time.time()
    print("Took {} seconds to execute.".format(b - a))
