import numpy as np
import pickle
import argparse
import rdkit.Chem.AllChem as rdkit
import itertools as it
import matplotlib.pyplot as plt


def get_fingerprint(molecule):
    return rdkit.GetMorganFingerprint(molecule.to_rdkit_mol(), 2)


def similarity(molecule1, molecule2):
    return rdkit.DataStructs.DiceSimilarity(
        bv1=get_fingerprint(molecule1),
        bv2=get_fingerprint(molecule2),
    )


def generation_similarity(generation):
    return np.mean([
        similarity(m1, m2)
        for m1, m2 in it.product(generation, generation)
    ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'progress_pkl',
        help='The path to the pickled EA progress.',
    )
    parser.add_argument(
        'output_file',
        help='The path to the plot.',
    )
    args = parser.parse_args()

    with open(args.progress_pkl, 'rb') as f:
        progress = pickle.load(f)

    fig, ax = plt.subplots()
    ys = list(map(generation_similarity, progress))
    xs = list(range(1, len(progress)+1))
    ax.plot(xs, ys)
    ax.set(xlabel='Generation', ylabel='Mean Dice Similarity')
    plt.savefig(args.output_file, dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
