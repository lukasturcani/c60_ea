import argparse
import os
from os.path import join
import pickle
import rdkit.Chem.AllChem as rdkit


def _write_fittest(progress, output_directory):
    for i, generation in enumerate(progress):
        rdkit.MolToMolFile(
            mol=generation[-1].cage,
            filename=join(
                output_directory,
                f'{i}_{generation[-1].name}.mol'
            ),
            kekulize=False,
            forceV3000=True,
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            'Write the MDL MOL file of the fittest molecule at each '
            'generation.'
        )
    )
    parser.add_argument(
        'progress_pkl',
        help='The path to the pickled EA progress.',
    )
    parser.add_argument(
        'output_directory',
        help='The path to the output directory.',
    )

    args = parser.parse_args()
    if os.path.exists(args.output_directory):
        os.rmdir(args.output_directory)
    os.makedirs(args.output_directory)

    with open(args.progress_pkl, 'rb') as f:
        progress = pickle.load(f)

    _write_fittest(progress, args.output_directory)


if __name__ == '__main__':
    main()
