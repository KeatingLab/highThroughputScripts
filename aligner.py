'''
Contains the Aligner class, which performs general alignment tasks such as
pairwise alignment, scoring, formatting, and iterating.
'''

NO_SCORE = -1e30

class Aligner(object):

    def __init__(self, mutation_threshold=-1, identical_score=1, different_score=-1, gap_score=0, gap_character=" "):
        self.mutation_threshold = mutation_threshold
        self.identical_score = identical_score
        self.different_score = different_score
        self.gap_score = gap_score
        self.gap_character = gap_character

    def scoring_map(self, length, ranges):
        '''
        From the given length and list of range tuples, returns a list of boolean values
        that indicate whether or not to use the base at each index in scoring.
        '''
        ret = [False for _ in xrange(length)]
        for start, end in ranges:
            for i in xrange(start, end):
                ret[i] = True
        return ret

    def score(self, sequence_1, sequence_2, offset, scoring_maps=None):
        '''
        Scores the alignment where sequence_2 is shifted by offset to the right
        of the start of sequence_1.
        '''
        total = 0
        num_mutations = 0
        if offset < 0:
            total += -offset * self.gap_score
        for i, base_1 in enumerate(sequence_1):
            if scoring_maps is not None and ((scoring_maps[0] is not None and not scoring_maps[0][i]) or (scoring_maps[1] is not None and not scoring_maps[1][i - offset])):
                total += self.gap_score
            elif i - offset < 0 or i - offset >= len(sequence_2):
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

    def length(self, *sequences):
        '''
        Takes tuples of the form (sequence, offset) and returns the total length
        of the alignment.
        '''
        if len(sequences) == 0: return 0
        min_offset = min(sequences, key=lambda x: x[1])[1]
        rescaled_sequences = [(sequence, offset - min_offset) for sequence, offset in sequences]

        return max(offset + len(sequence) for sequence, offset in rescaled_sequences)

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

    def align(self, sequence_1, sequence_2, min_overlap=-1, max_overlap=-1, unidirectional=False, reverse=False, scoring_ranges=None):
        '''
        If unidirectional is True, then sequence_2 will be forced to start at or
        after the start of sequence_1.

        If scoring_ranges is not None, it should be a tuple (ranges_1, ranges_2)
        where each element is a list of ranges of bases that should be scored, for
        sequence_1 and sequence_2, respectively.
        '''
        offset_range = xrange(0 if unidirectional else -len(sequence_2) + 1, len(sequence_1))
        if reverse:
            offset_range = reversed(offset_range)

        best_offset = 0
        best_score = NO_SCORE

        scoring_maps = [None, None]
        if scoring_ranges is not None:
            if scoring_ranges[0] is not None:
                scoring_maps[0] = self.scoring_map(len(sequence_1), scoring_ranges[0])
            if scoring_ranges[1] is not None:
                scoring_maps[1] = self.scoring_map(len(sequence_2), scoring_ranges[1])

        for offset in offset_range:
            overlap = self.overlap(sequence_1, sequence_2, offset)
            if (min_overlap >= 0 and overlap < min_overlap) or (max_overlap >= 0 and overlap > max_overlap):
                continue
            score = self.score(sequence_1, sequence_2, offset, scoring_maps=scoring_maps)
            if score != NO_SCORE and score > best_score:
                best_score = score
                best_offset = offset
        return (best_offset, best_score)

    def best_alignment(self, target, sequences, candidate_scoring_ranges=None, **kwargs):
        '''
        Returns (best_sequence_index, best_offset, best_score), where best_sequence_index
        is the index of the sequence in `sequences` that results in the best
        alignment with the target. If candidate_scoring_ranges is not None, it should be
        a list of scoring ranges for each sequence in `sequences`. Any other
        optional arguments are passed through to the `align` function.
        '''
        if candidate_scoring_ranges is None:
            candidate_scoring_ranges = [None for i in xrange(len(sequences))]
        results = [[i] + list(self.align(sequence, target, scoring_ranges=(candidate_scoring_ranges[i], None), **kwargs)) for i, sequence in enumerate(sequences)]
        return max(results, key=lambda x: x[2])
