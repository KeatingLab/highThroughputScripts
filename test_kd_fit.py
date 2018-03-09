'''
This script takes the inputs and outputs to TiteSeq as parameters, and provides
a series of statistical measures of the fit provided by the Kd fitting algorithm.
'''

import os
import numpy as np
import time
import stat_collector as sc
import argparse
import pandas as pd
from scipy.stats import pearsonr, chisquare

### TiteSeq constants
from titeseq_utils import *

### CSV keys

SEQUENCE_KEY = "sequence"
INPUT_X_CSV_PREFIX = "input_x_"
OUTPUT_X_CSV_PREFIX = "output_x_"
X_CORRELATION_R_KEY = "x correlation (R^2)"
X_CORRELATION_CHI_KEY = "x correlation (chi^2)"
X_AVG_BIN_CORRELATION_R_KEY = "x avg bin correlation (R^2)"
X_AVG_BIN_CORRELATION_CHI_KEY = "x avg bin correlation (chi^2)"
INPUT_X_KD_FIT_CORRELATION_KEY = "input x/Kd correlation"
OUTPUT_X_KD_FIT_CORRELATION_KEY = "output x/Kd correlation"

'''
These functions indicate how to obtain the fit x values, Kd, S, and b values for
the fit curve from a given row in the output CSV file.
'''
get_fit_x = lambda row: row.iloc[1:33]
get_Kd = lambda row: row.iloc[33]   # row.iloc[9] for naive
get_S = lambda row: row.iloc[34]    # row.iloc[10] for naive
get_B = lambda row: basal           # row.iloc[11] for naive

'''
List of expression strings which will be evaluated and written out to the
params.txt file.
'''
PARAMETER_LIST = ["args.ts_input", "args.ts_output", "args.output", "READ_COUNT_COLUMN_RANGE", "NORMALIZED_COUNT_COLUMN_RANGE", "BINS_VARY_FIRST", "BINS_ORDER_ASCENDING", "CONCENTRATIONS_ORDER_ASCENDING", "bin_fluorescences", "concentrations"]

def read_titeseq_items(input, output):
    '''
    Reads the given input and output files (to TiteSeq, not to this script), and
    returns a list of lists of information about each sequence. Each inner list
    will contain the following: the sequence, the original read counts array
    (as a numpy array), the read counts array predicted by the TiteSeq
    algorithm, and the parameters for the fitting curve produced by the
    algorithm (K, s, b, as a tuple).
    '''
    in_csv = pd.read_csv(input, header=None)
    out_csv = pd.read_csv(output, header=None)

    seq_info = {}
    for _, row in out_csv.iterrows():
        seq = row.iloc[0]
        seq_item = [seq]
        # Assume the TiteSeq output CSV is ordered correctly - ascending concentrations, then ascending bins within them
        fit_x = normalized_count_array(reshape_from_csv_format(get_fit_x(row).as_matrix().astype(np.float32), use_input_order=False))
        seq_item.append(fit_x)
        seq_item += [get_Kd(row), get_S(row), get_B(row)]
        seq_info[seq] = seq_item

    num_completed = 0
    raw_start, raw_end = READ_COUNT_COLUMN_RANGE
    col_start, col_end = NORMALIZED_COUNT_COLUMN_RANGE

    # Sum normalized counts per CSV column, to account for sequencing depth
    total_read_counts = in_csv.iloc[:,raw_start:raw_end].sum().as_matrix().astype(np.float32)
    T = reshape_from_csv_format(total_read_counts)

    for _, row in in_csv.iterrows():
        seq = row.iloc[0]
        if seq not in seq_info:
            continue
        seq_item = seq_info[seq]
        num_completed += 1
        original_counts = reshape_from_csv_format(row.iloc[col_start:col_end].as_matrix().astype(np.float32))
        seq_item.insert(1, normalized_count_array(original_counts / T))

    assert num_completed == len(seq_info), "Incomplete sequences found (sequences in TiteSeq output CSV but not in input)"

    return seq_info

def test_kd_fits(ts_input, ts_output, output):
    '''
    Reads the read count arrays from the given TiteSeq input/output files, and
    determines the fit between the actual and fit counts as well as the strength
    of correlation between the actual/predicted counts and the fit Kd curve.
    Writes the results to the given output path.
    '''
    seqs = read_titeseq_items(ts_input, ts_output)

    data = []
    for seq, info in seqs.iteritems():
        seq, in_R, out_R, K, S, b = info
        flat_in_R = in_R.flatten()
        flat_out_R = out_R.flatten()

        stats = {SEQUENCE_KEY: seq}
        for i in range(COUNT_ARRAY_SIZE):
            stats[INPUT_X_CSV_PREFIX + str(i)] = flat_in_R[i]
            stats[OUTPUT_X_CSV_PREFIX + str(i)] = flat_out_R[i]

        stats[X_CORRELATION_R_KEY] = pearsonr(flat_in_R, flat_out_R)[0]
        stats[X_CORRELATION_CHI_KEY] = chisquare(flat_in_R, flat_out_R)[0]

        avg_in_R = average_bin_positions(in_R)
        avg_out_R = average_bin_positions(out_R)
        stats[X_AVG_BIN_CORRELATION_R_KEY] = pearsonr(avg_in_R, avg_out_R)[0]
        stats[X_AVG_BIN_CORRELATION_CHI_KEY] = chisquare(avg_in_R, avg_out_R)[0]

        # Compute an "expected" distribution based on the Kd curve
        data.append(stats)


    cols = [SEQUENCE_KEY] + [INPUT_X_CSV_PREFIX + str(i) for i in range(COUNT_ARRAY_SIZE)] + [OUTPUT_X_CSV_PREFIX + str(i) for i in range(COUNT_ARRAY_SIZE)] + [X_CORRELATION_R_KEY, X_CORRELATION_CHI_KEY, X_AVG_BIN_CORRELATION_R_KEY, X_AVG_BIN_CORRELATION_CHI_KEY, INPUT_X_KD_FIT_CORRELATION_KEY, OUTPUT_X_KD_FIT_CORRELATION_KEY]

    if len(data) > 0:
        df = pd.DataFrame(data, columns=cols)
        df.to_csv(output, index=False, float_format="%.4e")

if __name__ == '__main__':
    a = time.time()  # Time the script started
    parser = argparse.ArgumentParser(description='Takes the inputs and outputs to TiteSeq as parameters, and provides a series of statistical measures of the fit provided by the Kd fitting algorithm.')
    parser.add_argument('ts_input', metavar='TI', type=str,
                        help='The path to the CSV input to TiteSeq')
    parser.add_argument('ts_output', metavar='TO', type=str,
                        help='The path to the CSV output to TiteSeq')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the directory to write output')
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    test_kd_fits(args.ts_input, args.ts_output, os.path.join(args.output, os.path.basename(args.ts_output)))

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))

    # Write out the input parameters to file
    params_path = os.path.join(args.output, "params.txt")
    sc.write_input_parameters([(name, eval(name)) for name in PARAMETER_LIST], params_path)
