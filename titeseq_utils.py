import numpy as np
from scipy.special import factorial

'''
The range of indexes of the columns in the TiteSeq input CSV to use as
normalized counts.
'''
NORMALIZED_COUNT_COLUMN_RANGE = (33, 65)

'''
The range of indexes of the columns in the TiteSeq input CSV to use as raw counts.
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
Parameters (init0, s0, Kd0) for Mass-Titer.
'''
MT_LOWER = 1.0
MT_UPPER = 4.0
MT_KD_INITIAL = 1.0e-7

'''
The minimum raw reads (over all columns) required to compute the KD of a sequence.
'''
MINIMUM_TOTAL_READ_COUNT = 100

'''
The number of total elements in any read count array (automatically calculated
based on concentrations and bin_fluorescences).
'''
COUNT_ARRAY_SIZE = concentrations.shape[0] * bin_fluorescences.shape[0]

def reshape_from_csv_format(matrix, use_input_order=True):
    '''
    Formats the given 1-D numpy array so that it fits the requirements for
    TiteSeq, given the user-provided settings above regarding the dimensions and
    ordering of the CSV values.

    If use_input_order is True, the matrix will be reshaped using the above
    settings. If it is False, the matrix will be assumed to be in the correct
    order already: c1b1, c1b2, c1b3, c1b4, c2b1, c2b2, ... where c is concentration
    and b is bin.
    '''
    # Get dimensions of matrices from user-provided arrays
    M = concentrations.shape[0]
    N = bin_fluorescences.shape[0]

    # Final dimensions should be M x N
    if use_input_order:
        if BINS_VARY_FIRST:
            matrix = matrix.reshape(M, N)
        else:
            matrix = matrix.reshape(N, M).T

        if not BINS_ORDER_ASCENDING:
            matrix = np.flip(matrix, 1)
        if not CONCENTRATIONS_ORDER_ASCENDING:
            matrix = np.flip(matrix, 0)
    else:
        matrix = matrix.reshape(M, N)

    return matrix

def create_b(R):
    '''
    Determines the b matrix from the R matrix.
    '''
    b = (R.sum(axis=1, keepdims=True) / R).reshape(8, 4)
    b[np.isnan(b) + np.isinf(b)] = DEFAULT_B_VALUE
    return b

def normalized_count_array(R):
    '''
    Normalizes the given count array so that each row sums to one (distribution
    over bin counts).
    '''
    return R / R.sum(axis=1, keepdims=True)

def average_bin_positions(R, already_normalized=False):
    '''
    Returns a vector containing the average bin positions for each concentration.
    '''
    if already_normalized:
        return R[:,0] * (len(bin_fluorescences) * (len(bin_fluorescences) + 1) / 2.0)
    return np.average(np.vstack([range(1, len(bin_fluorescences) + 1)] * len(concentrations)).astype(float), axis=1, weights=normalized_count_array(R))

def fake_bin_data_from_averages(averages):
    '''
    Returns a matrix that, when average_bin_positions is applied to it, will yield
    the correct averages.
    '''
    return np.tile(averages, (len(bin_fluorescences), 1)).T / (len(bin_fluorescences) * (len(bin_fluorescences) + 1) / 2.0)
