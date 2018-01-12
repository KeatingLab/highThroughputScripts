'''
The purpose of this script is to take a single file containing sequential pairs
of forward and reverse reads, and to align each pair together to produce a single
DNA sequence, protein sequence, and set of quality scores.
'''

from Bio import SeqIO,Seq
import multiprocessing
import sys, os
from itertools import izip
import math
from functools import partial
import time
import argparse

NO_SCORE = -1e30
FORMAT = 'fastq-sanger'
SEQUENCE_QUALITY_KEY = 'phred_quality'
UNSPECIFIED_BASE = 'N'

'''
These determine how statistics are collected when the --stats parameter is set.
The output will be in CSV format with three elements per row:

* The first number is a quality threshold (n / 10 gives the negative log probability
of a sequencing error).
* The second number is the number of bases that are allowed to fall below the
quality threshold.
* The third number is the number of reads that would pass the given threshold and
tolerance.
'''
STAT_CUTOFFS = [0, 1, 2, 3, 5, 10, 20]
STAT_QUALITY_BINS = xrange(0, 40, 5)

'''
The references are the sequences onto which the forward and reverse reads will be
scaffolded. The output will be indexed by the number of the reference to which
each read was aligned.
'''
REFERENCE_SEQUENCES = ["ACTCGTCCCAAACAAGAACCTCAGGAAATCGATTTCCCGGACGATCTGCCAAACGACTTATCACG"]

'''
If OUTPUT_RANGES is not None, it should be a list with the same length as
REFERENCE_SEQUENCES. When a pair of reads is aligned to one of the reference
sequences, the sequence that is output will be trimmed to the range given by this
setting.
'''
OUTPUT_RANGES = None

class Aligner(object):
    '''
    Performs general alignment tasks, such as pairwise alignment, scoring,
    formatting, and iterating.
    '''

    def __init__(self, mutation_threshold=-1, identical_score=1, different_score=-1, gap_score=0, gap_character=" "):
        self.mutation_threshold = mutation_threshold
        self.identical_score = identical_score
        self.different_score = different_score
        self.gap_score = gap_score
        self.gap_character = gap_character

    def score(self, sequence_1, sequence_2, offset):
        '''
        Scores the alignment where sequence_2 is shifted by offset to the right
        of the start of sequence_1.
        '''
        total = 0
        num_mutations = 0
        if offset < 0:
            total += -offset * self.gap_score
        for i, base_1 in enumerate(sequence_1):
            if i - offset < 0 or i - offset >= len(sequence_2):
                total += self.gap_score
            else:
                base_2 = sequence_2[i - offset]
                if base_1 == base_2:
                    total += self.identical_score
                elif self.mutation_threshold < 0 or num_mutations < self.mutation_threshold:
                    total += self.different_score
                    num_mutations += 1
                else:
                    return NO_SCORE
        if len(sequence_1) - offset < len(sequence_2):  # Still bases left in sequence_2
            total += (len(sequence_1) - offset) * self.gap_score
        return total

    def overlap(self, sequence_1, sequence_2, offset):
        '''
        Computes the number of bases that overlap between the two sequences when
        sequence_2 is shifted by offset to the right of the start of sequence_1.
        '''
        if offset < 0:
            # sequence_2 is hanging off the 5' end of sequence_1
            return min(len(sequence_1), len(sequence_2) + offset)
        elif offset + len(sequence_2) >= len(sequence_1):
            # sequence_2 is hanging off the 3' end of sequence_1
            return len(sequence_1) - offset
        else:
            # sequence_2 is contained within sequence_1
            return len(sequence_2)

    def format(self, sequence_1, sequence_2, offset):
        '''
        Returns a tuple (sequence_1, sequence_2) where each sequence string has
        been padded with spaces to be the same length, reflecting the given offset.
        '''
        return self.format_multiple((sequence_1, 0), (sequence_2, offset))

    def format_multiple(self, *sequences):
        '''
        Takes tuples of the form (sequence, offset) and returns a list of
        sequences where each sequence been padded with spaces to be the same
        length, reflecting the given offsets.
        '''
        if len(sequences) == 0: return []
        min_offset = min(sequences, key=lambda x: x[1])[1]
        rescaled_sequences = [(sequence, offset - min_offset) for sequence, offset in sequences]
        rets = [self.gap_character * offset + sequence for sequence, offset in rescaled_sequences]

        max_length = len(rets[0])
        while True:
            changed = False
            for i in xrange(len(rets)):
                if len(rets[i]) < max_length:
                    rets[i] += self.gap_character * (max_length - len(rets[i]))
                elif len(rets[i]) > max_length:
                    max_length = len(rets[i])
                    changed = True
            if not changed:
                break
        return rets

    def enumerate(self, sequence_1, sequence_2, offset, placeholder=None):
        '''
        Enumerates the alignment specified by offset and yields a tuple (base_1, base_2)
        for each position. One of the bases may be `placeholder` (or by default
        the gap character) if only one sequence is present at a position.

        Note that the sequences can be lists or strings; for instance, to enumerate
        the quality scores for a given alignment, simply pass the lists of quality
        scores for sequence_1 and sequence_2.
        '''
        for item in self.enumerate_multiple((sequence_1, 0), (sequence_2, offset), placeholder=placeholder):
            yield item

    def enumerate_multiple(self, *sequences):
        '''
        Enumerates the alignment specified by the given (sequence, offset) pairs and
        yields a list [base_1, base_2, ...] for each position. One of the bases may
        be None if only one sequence is present at a position.

        Note that the sequences can be lists or strings; for instance, to enumerate
        the quality scores for a given alignment, simply pass the lists of quality
        scores for sequence_1 and sequence_2.
        '''
        placeholder = None
        if len(sequences) == 0: return
        min_offset = min(sequences, key=lambda x: x[1])[1]
        rescaled_sequences = [(sequence, offset - min_offset) for sequence, offset in sequences]

        index = 0
        added_element = True
        while added_element:
            item = []
            added_element = False
            for sequence, offset in rescaled_sequences:
                if index - offset < len(sequence) and index - offset >= 0:
                    item.append(sequence[index - offset])
                    added_element = True
                else:
                    item.append(placeholder)
            index += 1
            if added_element:
                yield item

    def align(self, sequence_1, sequence_2, min_overlap=-1, max_overlap=-1, unidirectional=False, reverse=False):
        '''
        If unidirectional is True, then sequence_2 will be forced to start at or
        after the start of sequence_1.
        '''
        offset_range = xrange(0 if unidirectional else -len(sequence_2) + 1, len(sequence_1))
        if reverse:
            offset_range = reversed(offset_range)

        best_offset = 0
        best_score = NO_SCORE
        for offset in offset_range:
            overlap = self.overlap(sequence_1, sequence_2, offset)
            if (min_overlap >= 0 and overlap < min_overlap) or (max_overlap >= 0 and overlap > max_overlap):
                continue
            score = self.score(sequence_1, sequence_2, offset)
            if score != NO_SCORE and score > best_score:
                best_score = score
                best_offset = offset
        return (best_offset, best_score)

    def best_alignment(self, target, sequences, **kwargs):
        '''
        Returns (best_sequence_index, best_offset, best_score), where best_sequence_index
        is the index of the sequence in `sequences` that results in the best
        alignment with the target. Any optional arguments are passed through to
        the `align` function.
        '''
        results = [[i] + list(self.align(sequence, target, **kwargs)) for i, sequence in enumerate(sequences)]
        return max(results, key=lambda x: x[2])


def is_quality_read(record, threshold=20, allowable_misreads=0):
    '''
    Checks the quality of the given SeqRecord. If more than `allowable_misreads`
    bases have quality below `threshold`, this function returns False, else it
    returns True.
    '''
    misread_count = sum(1 for q in record.letter_annotations[SEQUENCE_QUALITY_KEY] if q < threshold)
    return misread_count <= allowable_misreads

def quality_stats(qualities, cutoffs=[0, 2, 4, 6, 10]):
    '''
    Determines and returns the minimum score when the given cutoff numbers of bad
    reads are excluded.
    '''
    sorted_qualities = sorted(qualities)
    return list(map(lambda x: sorted_qualities[x], cutoffs))

def combine_records(forward_record, reverse_record, reference_sequences, min_overlap=-1, max_overlap=-1):
    '''
    Computes the alignments of both forward and reverse reads to the reference
    sequences. Synthesizes those alignments, using the better-quality read in
    the case of a conflict. Returns (index, sequence, quality) where `index` is
    the index of the reference sequence used, `sequence` is the combined DNA
    sequence, and `quality` is the quality of each base in the combined sequence.

    The optional parameters min_overlap and max_overlap correspond to the overlap
    constraints on the alignment between the forward and reverse reads.
    '''
    aligner = Aligner()

    forward_str = str(forward_record.seq)
    reverse_str = str(reverse_record.seq.reverse_complement())

    reference_index, forward_offset, forward_score = aligner.best_alignment(forward_str, reference_sequences, unidirectional=True, min_overlap=len(forward_str))
    reverse_offset, _ = aligner.align(forward_str, reverse_str, unidirectional=True, reverse=True, min_overlap=10, max_overlap=30)
    reverse_offset += forward_offset # Accounts for alignment to reference
    reference = reference_sequences[reference_index]

    combined_sequence = ""
    combined_quality = []

    # The aligner will enumerate the aligned characters or elements of each iterable we give it.
    # Zipping generators for both the sequence and the quality allows us to enumerate them together.
    sequence_generator = aligner.enumerate_multiple((reference, 0),
                                                    (forward_str, forward_offset),
                                                    (reverse_str, reverse_offset))
    quality_generator = aligner.enumerate_multiple(([None for i in xrange(len(reference))], 0),
                                                   (forward_record.letter_annotations[SEQUENCE_QUALITY_KEY], forward_offset),
                                                   (reverse_record.letter_annotations[SEQUENCE_QUALITY_KEY], reverse_offset))
    for bases, qualities in izip(sequence_generator, quality_generator):
        _, forward_base, reverse_base = bases
        _, forward_quality, reverse_quality = qualities

        if forward_base is None and reverse_base is None:
            combined_sequence += UNSPECIFIED_BASE
            combined_quality.append(0)
        elif forward_base is None:
            combined_sequence += reverse_base
            combined_quality.append(reverse_quality)
        elif reverse_base is None:
            combined_sequence += forward_base
            combined_quality.append(forward_quality)
        else:
            base, quality = max([(forward_base, forward_quality), (reverse_base, reverse_quality)], key=lambda x: x[1])
            combined_sequence += base
            combined_quality.append(quality)

    return reference_index, combined_sequence, combined_quality

def combine_records_processor(references, (forward, reverse), threshold=0, output_ranges=None, **kwargs):
    '''
    Performs initial quality check and passes through to the main combine_records
    function. Returns (input, reference, dna_sequence, aa_sequence, quality), where reference
    is the index of the used reference. If an alignment was not computed because of
    bad quality, reference will be -1.

    If output_ranges is not None, its element corresponding to the chosen reference
    sequence will define the range of bases in the scaffolding sequence that
    should be output. (Note: it should match the reading frame of the coding
    sequence to ensure correct translation.)
    '''
    if not is_quality_read(forward, threshold=threshold) or not is_quality_read(reverse, threshold=threshold):
        return -1, None, None, None
    reference, dna_sequence, quality = combine_records(forward, reverse, references, **kwargs)
    if len(dna_sequence) == 0:
        return -1, None, None, None

    if output_ranges is not None:
        start, end = output_ranges[reference]
        dna_sequence = dna_sequence[start:end]
        quality = quality[start:end]

    aa_sequence = str(Seq.Seq(dna_sequence).translate())
    return (forward, reverse), reference, dna_sequence, aa_sequence, quality

def update_quality_stats(quality_dict, record):
    '''
    Helper function for write_combined_records that updates the given quality
    dictionary with qualities provided in the given SeqRecord (or list) object.
    '''
    quality_list = record.letter_annotations[SEQUENCE_QUALITY_KEY] if type(record) == SeqIO.SeqRecord else record
    stats = quality_stats(quality_list, STAT_CUTOFFS)
    for cutoff, stat in izip(STAT_CUTOFFS, stats):
        for quality in STAT_QUALITY_BINS:
            if quality > stat: break
            quality_dict[quality][cutoff] += 1

def write_combined_records(input_path, references, out_dir, num_processes=15, threshold=0, stats=False, output_ranges=None):
    '''
    Combines the records at input_path by aligning them to the given reference sequences,
    and saves them to the appropriate locations within out_dir. The combined DNA
    sequences are written to out_dir/dnaframe, the amino acid sequences are
    written to out_dir/seqframe, and the qualities are written to qualframe (as CSV).
    '''
    basename = os.path.basename(input_path)

    # Create directories if necessary
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    dna_path = os.path.join(out_dir, "dnaframe")
    aa_path = os.path.join(out_dir, "seqframe")
    qual_path = os.path.join(out_dir, "qualframe")
    for path in [dna_path, aa_path, qual_path]:
        if not os.path.exists(path):
            os.mkdir(path)
    if stats:
        stats_path = os.path.join(out_dir, "stats")
        if not os.path.exists(stats_path):
            os.mkdir(stats_path)

    # Open file streams
    dna_files = [open(os.path.join(dna_path, basename + "_{}".format(ref)), "w") for ref in xrange(len(references))]
    aa_files = [open(os.path.join(aa_path, basename + "_{}".format(ref)), "w") for ref in xrange(len(references))]
    qual_files = [open(os.path.join(qual_path, basename + "_{}".format(ref)), "w") for ref in xrange(len(references))]

    if stats:
        # Initialize statistics bookkeeping
        forward_quality_stats = {qual: {cutoff: 0 for cutoff in STAT_CUTOFFS} for qual in STAT_QUALITY_BINS}
        reverse_quality_stats = {qual: {cutoff: 0 for cutoff in STAT_CUTOFFS} for qual in STAT_QUALITY_BINS}
        total_quality_stats = {qual: {cutoff: 0 for cutoff in STAT_CUTOFFS} for qual in STAT_QUALITY_BINS}

    with open(input_path, 'rU') as file:

        records = SeqIO.parse(file, FORMAT)

        pool = multiprocessing.Pool(processes=num_processes)
        processor = partial(combine_records_processor,
                            references,
                            threshold=threshold,
                            output_ranges=output_ranges,
                            min_overlap=10,
                            max_overlap=30)
        for result in pool.imap(processor, izip(records, records), chunksize=1000):
            original_input, ref_index, dna_sequence, aa_sequence, quality = result
            if ref_index == -1:
                continue

            dna_files[ref_index].write(dna_sequence + "\n")
            aa_files[ref_index].write(aa_sequence + "\n")
            qual_files[ref_index].write(",".join([str(q) for q in quality]) + "\n")

            if stats:
                forward, reverse = original_input
                update_quality_stats(forward_quality_stats, forward)
                update_quality_stats(reverse_quality_stats, reverse)
                update_quality_stats(total_quality_stats, quality)

    for file in dna_files + aa_files + qual_files:
        file.close()
    if stats:
        write_quality_stats(forward_quality_stats, os.path.join(out_dir, "stats", basename + "_forward.txt"))
        write_quality_stats(reverse_quality_stats, os.path.join(out_dir, "stats", basename + "_reverse.txt"))
        write_quality_stats(total_quality_stats, os.path.join(out_dir, "stats", basename + "_total.txt"))

def write_quality_stats(quality_stats, out_path):
    with open(out_path, "w") as stats_file:
        for cutoff in sorted(quality_stats.keys()):
            for quality in sorted(quality_stats[cutoff].keys()):
                stats_file.write(",".join([str(cutoff), str(quality), str(quality_stats[cutoff][quality])]) + "\n")

if __name__ == '__main__':
    a = time.time()  # Time the script started
    parser = argparse.ArgumentParser(description='Reads a single barcode file containing pairs of forward and reverse reads, and aligns them to a set of possible scaffolds (defined at the top of this script). Writes the results to three directories, containing the DNA sequence, the AA sequence, and the quality scores, respectively.')
    parser.add_argument('input', metavar='I', type=str,
                        help='The path to the FASTQ input file')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the output directory')
    parser.add_argument('-t', '--threshold', type=int, default='0',
                        help='The quality score below which reads should not be used (default 0)')
    parser.add_argument('-p', '--processes', type=int, default=15,
                        help='The number of processes to use')
    parser.add_argument('--stats', dest='stats', action='store_true',
                        help='Whether to output stats into output/[base]_stats.txt')
    parser.set_defaults(stats=False)
    args = parser.parse_args()

    write_combined_records(args.input, REFERENCE_SEQUENCES, args.output, num_processes=args.processes, threshold=args.threshold, stats=args.stats, output_ranges=OUTPUT_RANGES)

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))
    #offset, score = aligner.align(seq2, seq1, unidirectional=True)
    #print('\n'.join(aligner.format(seq2, seq1, offset)))