import argparse
import os
import sklearn
import numpy as np
from permute_sequences import amino_acids
from sklearn.linear_model import LinearRegression
import sys

def one_hot(sequence):
    """
    Default feature transformation; takes the given sequence and encodes it into a feature vector by one-hot encoding.
    """
    def aa_vector(letter):
        vector = np.zeros(len(amino_acids))
        vector[amino_acids.index(letter)] = 1
        return vector
    return np.hstack([aa_vector(c) for c in sequence])

def score_sequences(model, sequences, outpath=None, feature_transform=one_hot, filter=None):
    """Computes the model's prediction for all sequences in the given iterable, and writes the (optionally filtered) results to outpath.

    Parameters:
    model -- A scikit-learn model with a predict function that can be used to compute a score or classification for the given sequences (transformed into feature vectors).
    sequences -- An iterable or generator that provides the peptide sequences for which to predict scores.
    outpath -- The destination file path for writing the results. The output will be in CSV format, where each row contains two elements: the sequence and the predicted score or classification. If outpath is None, writes to stdout.
    feature_transform -- A function that takes in a sequence and returns a numpy array representing the feature vector for that sequence. The default is to use one-hot encoding over the possible amino acids in each position.
    filter -- A function that takes a sequence and score as input and returns a Boolean value indicating whether or not to write the result to file. If this function is None, all results are written to file.
    """
    if outpath is not None:
        file = open(outpath, 'w')
    else:
        file = sys.stdout

    for seq in sequences:
        score = model.predict(feature_transform(seq))
        if filter is not None and not filter(seq, score):
            continue
        file.write("{},{}\n".format(seq, score))

    if outpath is not None:
        file.close()

if __name__ == "__main__":
    # Test data
    X = np.random.rand(100, 1)
    Y = 50.0 * X + 5.0 * (np.random.rand(100, 1) - 0.5)
    model = LinearRegression()
    model.fit(X, Y)
    print("All:")
    score_sequences(model, [3, 6, 9, 12, 15], feature_transform=lambda x: x)
    print("Filtered:")
    score_sequences(model, [3, 6, 9, 12, 15], feature_transform=lambda x: x, filter=lambda s, v: v > 200)
