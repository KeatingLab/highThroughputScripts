import argparse
import random
import os

amino_acids = "ARNDCEQGHILKMFPSTWYV"
bases = "AGCT"

def detect_language(template):
    for char in template:
        if char not in amino_acids and char in bases:
            print("Inferring DNA.")
            return bases
        elif char not in bases and char in amino_acids:
            print("Inferring protein.")
            return amino_acids
    assert False, "Can't detect whether to use amino acids or nucleotides"

def write_permutations_recursive(template, base, language, outfiles):
    if len(template) == 0:
        outfiles[0].write(base + '\n')
    elif template[0] == '*':
        for i, char in enumerate(language):
            if len(base) == 0:
                print("{} / {} first positions".format(i, len(language)))
            if len(base) == 1:
                outfile = outfiles[i % len(outfiles)]
                write_permutations_recursive(template[1:], base + char, language, [outfile])
            else:
                write_permutations_recursive(template[1:], base + char, language, outfiles)
    else:
        write_permutations_recursive(template[1:], base + template[0], language, outfiles)

def write_permutations(template, out_path, mode=None, split=1):
    if mode == 'aa':
        language = amino_acids
    elif mode == 'dna':
        language = bases
    else:
        language = detect_language(template)

    files = []
    if split > 1:
        base, ext = os.path.splitext(out_path)
        for i in range(split):
            files.append(open(base + '_{}'.format(i) + ext, 'w'))
    else:
        files.append(open(out_path, 'w'))

    write_permutations_recursive(template, "", language, files)

    for file in files:
        file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Writes all permutations of either DNA or protein sequences with the given template to file.')
    parser.add_argument('seq', metavar='S', type=str,
                        help='The template sequence, where * denotes a variable position')
    parser.add_argument('output', metavar='O', type=str,
                        help='The path to the output file')
    parser.add_argument('-m', '--mode', type=str, default=None,
                        help='The set of characters to permute (aa or dna), optional')
    parser.add_argument('-s', '--split', type=int, default=1,
                        help='The number of files to split the output into')
    args = parser.parse_args()

    write_permutations(args.seq, args.output, args.mode, args.split)
