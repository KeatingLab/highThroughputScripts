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
from run_titeseq import MODE_X_STAR, MODE_MT

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
INPUT_KD_KEY = "input Kd"
INPUT_S_KEY = "input S"
INPUT_B_KEY = "input B"
OUTPUT_KD_KEY = "output Kd"
OUTPUT_S_KEY = "output S"
OUTPUT_B_KEY = "output B"
INPUT_ERROR_KEY = "input error metric"
OUTPUT_ERROR_KEY = "output error metric"

COLUMN_ORDER = [SEQUENCE_KEY] + [INPUT_X_CSV_PREFIX + str(i) for i in range(COUNT_ARRAY_SIZE)] + [OUTPUT_X_CSV_PREFIX + str(i) for i in range(COUNT_ARRAY_SIZE)] + [X_CORRELATION_R_KEY, X_CORRELATION_CHI_KEY, X_AVG_BIN_CORRELATION_R_KEY, X_AVG_BIN_CORRELATION_CHI_KEY, INPUT_X_KD_FIT_CORRELATION_KEY, OUTPUT_X_KD_FIT_CORRELATION_KEY] + [OUTPUT_KD_KEY, OUTPUT_S_KEY, OUTPUT_B_KEY, OUTPUT_ERROR_KEY, INPUT_KD_KEY, INPUT_S_KEY, INPUT_B_KEY, INPUT_ERROR_KEY]

'''
These functions indicate how to obtain the fit x values, Kd, S, and b values for
the fit curve from a given row in the output CSV file.
'''
get_fit_x = lambda row: row.iloc[1:33]
get_Kd = lambda row: row.iloc[33]
get_S = lambda row: row.iloc[34]
get_B = lambda row: row.iloc[35]

x_star_get_kd_sigma = lambda row: row.iloc[36]
mt_get_r2 = lambda row: row.iloc[38]

'''
List of expression strings which will be evaluated and written out to the
params.txt file.
'''
PARAMETER_LIST = ["args.ts_input", "args.ts_output", "args.output", "READ_COUNT_COLUMN_RANGE", "NORMALIZED_COUNT_COLUMN_RANGE", "BINS_VARY_FIRST", "BINS_ORDER_ASCENDING", "CONCENTRATIONS_ORDER_ASCENDING", "bin_fluorescences", "concentrations", "args.both_ts_output"]

def fit_parameters_from_row(row, method):
    '''
    Extracts the fit parameters from the given row, including additional error
    values depending on the method used to generate this data.
    '''
    base = [get_Kd(row), get_S(row), get_B(row)]
    if method == MODE_X_STAR:
        base.append(x_star_get_kd_sigma(row))
    elif method == MODE_MT:
        base.append(mt_get_r2(row))
    return base

def read_titeseq_items(input, output, both_ts_output=False, input_method=None, output_method=None):
    '''
    Reads the given input and output files (to TiteSeq, not to this script), and
    returns a list of lists of information about each sequence. Each inner list
    will contain the following: the sequence, the original read counts array
    (as a numpy array), the read counts array predicted by the TiteSeq
    algorithm, and the parameters for the fitting curve produced by the
    algorithm (K, s, b, and any relevant error values). If both_ts_output is True,
    the return value contains an additional element, which is the (K, s, b) fit
    values from the ts_input data.
    '''
    in_csv = pd.read_csv(input, header=None)
    out_csv = pd.read_csv(output, header=None)

    seq_info = {}
    for _, row in out_csv.iterrows():
        seq = row.iloc[0]
        seq_item = [seq]
        # Assume the TiteSeq output CSV is ordered correctly - ascending concentrations, then ascending bins within them
        fit_x = reshape_from_csv_format(get_fit_x(row).as_matrix().astype(np.float32), use_input_order=False)
        seq_item.append(fit_x)
        seq_item.append(fit_parameters_from_row(row, output_method))
        seq_info[seq] = seq_item

    num_completed = 0
    raw_start, raw_end = READ_COUNT_COLUMN_RANGE
    if both_ts_output:
        # The ts_input file is a run_titeseq output file. Get the x array from positions 1-NM.
        raw_start = 1
        raw_end = 1 + raw_end - raw_start
        col_start, col_end = raw_start, raw_end
        T = reshape_from_csv_format(np.ones(raw_end - raw_start))
    else:
        # The ts_input file is a run_titeseq input file, which has both raw and normalized counts
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
        original_counts = reshape_from_csv_format(row.iloc[col_start:col_end].as_matrix().astype(np.float32), use_input_order=not both_ts_output)
        seq_item.insert(1, original_counts / T)

        if both_ts_output:
            seq_item.append(fit_parameters_from_row(row, input_method))

    assert num_completed == len(seq_info), "Incomplete sequences found (sequences in TiteSeq output CSV but not in input)"

    return seq_info

def test_kd_fits(ts_input, ts_output, output, both_ts_output=False, input_method=None, output_method=None):
    '''
    Reads the read count arrays from the given TiteSeq input/output files, and
    determines the fit between the actual and fit counts as well as the strength
    of correlation between the actual/predicted counts and the fit Kd curve.
    Writes the results to the given output path.
    '''
    seqs = read_titeseq_items(ts_input, ts_output, both_ts_output, input_method, output_method)

    data = []
    for seq, info in seqs.iteritems():
        if both_ts_output:
            seq, in_R, out_R, fit, in_fit = info
        else:
            seq, in_R, out_R, fit = info

        flat_in_R = in_R.flatten()
        flat_out_R = out_R.flatten()

        stats = {SEQUENCE_KEY: seq}
        for i in range(COUNT_ARRAY_SIZE):
            stats[INPUT_X_CSV_PREFIX + str(i)] = flat_in_R[i]
            stats[OUTPUT_X_CSV_PREFIX + str(i)] = flat_out_R[i]

        avg_in_R = average_bin_positions(in_R, already_normalized=(both_ts_output and input_method == MODE_MT))
        avg_out_R = average_bin_positions(out_R, already_normalized=(output_method == MODE_MT))
        stats[X_AVG_BIN_CORRELATION_R_KEY] = pearsonr(avg_in_R, avg_out_R)[0] ** 2
        stats[X_AVG_BIN_CORRELATION_CHI_KEY] = chisquare(avg_in_R, avg_out_R)[0]

        if both_ts_output:
            stats[INPUT_KD_KEY] = in_fit[0]
            stats[INPUT_S_KEY] = in_fit[1]
            stats[INPUT_B_KEY] = in_fit[2]
            if len(in_fit) > 3:
                stats[INPUT_ERROR_KEY] = in_fit[3]
        stats[OUTPUT_KD_KEY] = fit[0]
        stats[OUTPUT_S_KEY] = fit[1]
        stats[OUTPUT_B_KEY] = fit[2]
        if len(fit) > 3:
            stats[OUTPUT_ERROR_KEY] = fit[3]
        # TODO: Compute an "expected" distribution based on the Kd curve

        flat_in_R = normalized_count_array(in_R).flatten()
        flat_out_R = normalized_count_array(out_R).flatten()
        stats[X_CORRELATION_R_KEY] = pearsonr(flat_in_R, flat_out_R)[0] ** 2
        stats[X_CORRELATION_CHI_KEY] = chisquare(flat_in_R, flat_out_R)[0]

        data.append(stats)

    if len(data) > 0:
        df = pd.DataFrame(data, columns=COLUMN_ORDER)
        df.to_csv(output, index=False, float_format="%.4e")

if __name__ == '__main__':
    a = time.time()  # Time the script started
    parser = argparse.ArgumentParser(description='Takes the inputs and outputs to TiteSeq as parameters, and provides a series of statistical measures of the fit provided by the Kd fitting algorithm.')
    parser.add_argument('ts_input', metavar='TI', type=str,
                        help='The path to the CSV input to TiteSeq')
    parser.add_argument('ts_output', metavar='TO', type=str,
                        help='The path to the CSV output to TiteSeq')
    parser.add_argument('-im', '--input-method', type=str, default=None,
                        dest='input_method', help='The mode in which run_titeseq was run to create the ts_input')
    parser.add_argument('-om', '--output-method', type=str, default=None,
                        dest='output_method', help='The mode in which run_titeseq was run to create the ts_output')
    parser.add_argument('-bto', '--both-ts-output', action='store_true',
                        dest='both_ts_output', help='Provide this flag if the ts_input file is formatted like an output file from run_titeseq.')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the CSV file in which to write output')
    args = parser.parse_args()

    if not os.path.exists(os.path.dirname(args.output)):
        os.mkdir(os.path.dirname(args.output))

    test_kd_fits(args.ts_input, args.ts_output, args.output, both_ts_output=args.both_ts_output, input_method=args.input_method, output_method=args.output_method)

    b = time.time()
    print("Took {} seconds to execute.".format(b - a))

    # Write out the input parameters to file
    params_path = os.path.join(os.path.dirname(args.output), "params.txt")
    sc.write_input_parameters([(name, eval(name)) for name in PARAMETER_LIST], params_path)
