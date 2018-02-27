'''
This script takes the CSV input generated by 05_write_titeseq_input and performs
the Tite-Seq analysis in parallel.
'''

from KD_fit_log_poiss import x_star
import numpy as np
import time
import stat_collector as sc
import argparse
import os
import pandas as pd
from itertools import imap
import multiprocessing

'''
The range of indexes of the columns in the CSV to use as normalized counts.
'''
NORMALIZED_COUNT_COLUMN_RANGE = (33, 65)

'''
The range of indexes of the columns in the CSV to use as raw counts.
'''
READ_COUNT_COLUMN_RANGE = (1, 33)

'''
Indicates whether the data groups values by concentrations first or by bins first.
'''
BINS_VARY_FIRST = True

'''
Indicates whether or not to flip the CSV values to match the order of bin_fluorescences
and concentrations below.
'''
BINS_ORDER_ASCENDING = False
CONCENTRATIONS_ORDER_ASCENDING = False

'''
Defines the values for the mean bin fluorescences (c) and the concentrations (fl).
The lengths of these arrays defines the dimensions of the matrix that will be used
for the input to TiteSeq. These values should be provided in *ascending* order,
regardless of whether or not the CSV counts are flipped.
'''
bin_fluorescences = np.array([3.0, 250.0, 900.0, 3000.0])
concentrations = np.array([2e-9, 6e-9, 2e-8, 6e-8, 2e-7, 6e-7, 2e-6, 6e-6])

'''
Miscellaneous other parameters for x_star (see KD_fit_log_poiss.py)
'''
k_scan = np.array([5e-5, 5e-4, 5e-3])
basal = 0.1
KD_scan = np.array(np.logspace(-5, -10, 71))
s_scan = (np.logspace(2, np.log10(bin_fluorescences[-1] - basal) - 0.02, 70))

'''
The default value to use in the b matrix if the value obtained is NaN or inf.
'''
DEFAULT_B_VALUE = 10.0

'''
The minimum raw reads (over all columns) required to compute the KD of a sequence.
'''
MINIMUM_TOTAL_READ_COUNT = 100

'''
List of expression strings which will be evaluated and written out to the
params.txt file.
'''
PARAMETER_LIST = ["args.input", "args.output", "READ_COUNT_COLUMN_RANGE", "NORMALIZED_COUNT_COLUMN_RANGE", "BINS_VARY_FIRST", "BINS_ORDER_ASCENDING", "CONCENTRATIONS_ORDER_ASCENDING", "bin_fluorescences", "concentrations", "k_scan", "basal", "KD_scan", "s_scan", "DEFAULT_B_VALUE", "MINIMUM_TOTAL_READ_COUNT"]

def reshape_from_csv_format(matrix):
    '''
    Formats the given 1-D numpy array so that it fits the requirements for
    TiteSeq, given the user-provided settings above regarding the dimensions and
    ordering of the CSV values.
    '''
    # Get dimensions of matrices from user-provided arrays
    M = concentrations.shape[0]
    N = bin_fluorescences.shape[0]

    # Final dimensions should be M x N
    if BINS_VARY_FIRST:
        matrix = matrix.reshape(M, N)
    else:
        matrix = matrix.reshape(N, M).T

    if not BINS_ORDER_ASCENDING:
        matrix = np.flip(matrix, 1)
    if not CONCENTRATIONS_ORDER_ASCENDING:
        matrix = np.flip(matrix, 0)

    return matrix

def create_b(R):
    '''
    Determines the b matrix from the R matrix.
    '''
    b = (R.sum(axis=1, keepdims=True) / R).reshape(8, 4)
    b[np.isnan(b) + np.isinf(b)] = DEFAULT_B_VALUE
    return b

def format_for_csv(value):
    '''Returns a string (in scientific notation) for the given number.'''
    return "{:.4e}".format(value)

def format_log_for_csv(value):
    '''Returns a shorter decimal string for the given number.'''
    return "{:.2f}".format(value)

def run_x_star_processor(input_vals):
    '''
    Helper function for run_x_star that runs in a multiprocessing pool. Takes
    an input of the form (R, T, b, seq) and returns (seq, result) where result
    is a tuple of the return values from x_star.
    '''
    if input_vals is None:
        return None
    R, T, b, seq = input_vals
    return seq, x_star(R, T, b, k_scan, bin_fluorescences, basal, KD_scan, s_scan, concentrations)

def run_x_star(in_path, out_path, line_num=None, num_processes=5):
    '''
    Takes as input the path to the output of 05_write_titeseq_input, which is a
    CSV containing sequences, their counts in all the Tite-Seq bins, and their
    normalized counts. This method interprets the normalized counts into the
    correct inputs for x_star (see KD_fit_log_poiss.py), runs x_star, and writes
    the output in CSV format to out_path.

    If line_num is not None, it can be an index representing the line number of the
    sequence for which to compute x_star.
    '''
    df = pd.read_csv(in_path, header=None)
    df = df.rename(columns={0:'seq'})

    # Clear output file
    if out_path is not None:
        open(out_path, 'w').close()

    raw_start, raw_end = READ_COUNT_COLUMN_RANGE
    col_start, col_end = NORMALIZED_COUNT_COLUMN_RANGE

    # Generate T matrix (sum of normalized counts per CSV column)
    total_read_counts = df.iloc[:,col_start:col_end].sum().as_matrix().astype(np.float32)
    T = reshape_from_csv_format(total_read_counts)

    def input_matrices(row_info):
        '''Generate R and b matrices from the given index and row.'''
        index, row = row_info
        if row.as_matrix(columns=df.columns[raw_start:raw_end]).sum() < MINIMUM_TOTAL_READ_COUNT:
            return None
        print("Processing sequence {} ({})...".format(index, row['seq']))
        R = reshape_from_csv_format(row.as_matrix(columns=df.columns[col_start:col_end]).astype(np.float32))
        b = create_b(R)
        return R, T, b, row['seq']

    def write_x_star_results(seq, result):
        '''Write the results of x_star to the output path or to the console.'''
        x, KD, s, KD_sigma, get_obj, prob, log_prob, try_LL, k_guess = result

        if out_path is not None:
            # Write results to CSV
            with open(out_path, 'a') as out_file:
                info_to_write = [seq] + [format_for_csv(val) for val in x.flatten()] + [format_for_csv(KD), format_for_csv(s), format_for_csv(KD_sigma)] + [format_for_csv(try_LL), format_for_csv(k_guess)]
                out_file.write(','.join(info_to_write) + '\n')
        else:
            # Print results to console
            print("{}: KD = {}".format(seq, KD))
            print("x = {}".format(x))

    if line_num is not None:
        # Process only the given line number
        row = df.iloc[line_num,:]
        print("Processing sequence {} ({})...".format(line_num, row[0]))
        write_x_star_results(*run_x_star_processor(input_matrices(row)))
    else:
        # Process all lines using a multiprocess pool
        pool = multiprocessing.Pool(processes=num_processes)
        inputs = imap(input_matrices, df.iterrows())
        for result in pool.imap(run_x_star_processor, inputs):
            if result is None:
                continue
            write_x_star_results(*result)


if __name__ == '__main__':
    a = time.time()  # Time the script started
    parser = argparse.ArgumentParser(description='Runs the TiteSeq calculation procedure on each of the sequences whose counts are given in the input file.')
    parser.add_argument('input', metavar='I', type=str,
                        help='The path to the input sequence counts directory')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the output directory')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    run_x_star(args.input, os.path.join(args.output, os.path.basename(args.input)))

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))

    # Write out the input parameters to file
    params_path = os.path.join(args.output, "params.txt")
    sc.write_input_parameters([(name, eval(name)) for name in PARAMETER_LIST], params_path)
