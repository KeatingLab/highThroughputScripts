'''
Reads a pair of forward and reverse files and extracts random entries according to
a specified probability. Run with Python 3.

Execution:
$ python random_sample_data.py forward_path reverse_path output_path -p 0.01
...
Generated test data at forward.fastq and reverse.fastq with <number> reads each
'''
import argparse
import random
import os, sys

# Number of lines per read
READ_LINE_COUNT = 4

PRINT_INTERVAL = 500000

def write_random_lines(forward_in, reverse_in, forward_out, reverse_out, proportion=0.005):
    """Reads lines from the input file handles and writes random groups of
    READ_LINE_COUNT lines to the output file handles. Returns the number of
    reads written."""

    i = 0
    should_write = False
    num_written = 0
    for f_line, r_line in zip(forward_in, reverse_in):
        if i % PRINT_INTERVAL == 0:
            print("Line {}, collected {} reads so far...".format(i, num_written))
        if i % 4 == 0:
            should_write = random.random() <= proportion
            if should_write:
                num_written += 1

        if should_write:
            forward_out.write(f_line)
            reverse_out.write(r_line)
        i += 1
    return num_written

if __name__ == '__main__':
    assert sys.version_info >= (3,0)

    parser = argparse.ArgumentParser(description='Shrink the ')
    parser.add_argument('forward', metavar='F', type=str,
                        help='The path to the forward reads')
    parser.add_argument('reverse', metavar='R', type=str,
                        help='The path to the reverse reads')
    parser.add_argument('out', metavar='O', type=str,
                        help='The path to the output directory')
    parser.add_argument('-p', '--proportion', type=float, default=0.005,
                        help='The approximate proportion of reads to include')
    args = parser.parse_args()

    f_file = open(args.forward, 'r')
    r_file = open(args.reverse, 'r')
    f_out = open(os.path.join(args.out, 'forward.fastq'), 'w')
    r_out = open(os.path.join(args.out, 'reverse.fastq'), 'w')

    num_written = write_random_lines(f_file, r_file, f_out, r_out, proportion=args.proportion)

    f_file.close()
    r_file.close()
    f_out.close()
    r_out.close()

    print("Generated test data at forward.fastq and reverse.fastq with {} reads each".format(num_written))
