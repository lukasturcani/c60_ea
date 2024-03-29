import os
import numpy as np
import pickle
import argparse
import rdkit.Chem.AllChem as rdkit
import itertools as it
import matplotlib.pyplot as plt


def get_fingerprint(molecule):
    return rdkit.GetMorganFingerprint(molecule, 2)


def similarity(molecule1, molecule2):
    return rdkit.DataStructs.DiceSimilarity(
        get_fingerprint(molecule1),
        get_fingerprint(molecule2),
    )


def get_building_blocks(generation):
    for cage in generation:
        for building_block in cage.building_blocks:
            rdkit.SanitizeMol(building_block.mol)
            yield building_block.mol


def dedupe(iterable):
    seen = set()
    for molecule in iterable:
        key = rdkit.MolToInchi(molecule)
        if key not in seen:
            seen.add(key)
            yield molecule


def generation_similarity(generation):
    building_blocks = list(dedupe(get_building_blocks(generation)))
    return np.mean([
        similarity(m1, m2)
        for m1, m2 in it.product(building_blocks, building_blocks)
    ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'progress_pkl',
        help='The path to the pickled EA progress.',
        nargs='+',
    )
    parser.add_argument(
        'output_file',
        help='The path to the plot.',
    )
    args = parser.parse_args()
    fig, ax = plt.subplots()

    for progress_pkl in args.progress_pkl:
        run, _ = os.path.splitext(os.path.basename(progress_pkl))
        with open(progress_pkl, 'rb') as f:
            progress = pickle.load(f)

        ys = list(map(generation_similarity, progress))
        xs = list(range(1, len(progress)+1))
        ax.plot(xs, ys, label=f'Run {run}')
        ax.set(
            xlabel='Generation',
            ylabel='Mean Dice Similarity',
        )

    ax.legend()
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(args.output_file, dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
