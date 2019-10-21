import itertools as it
import argparse
import stk
import numpy as np
import rdkit.Chem.AllChem as rdkit
import os
from os.path import join
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle


def _get_percent_rotatable_bonds(bb):
    return (
        100 *
        rdkit.CalcNumRotatableBonds(bb.mol) / bb.mol.GetNumBonds()
    )


def _get_percent_double_bonds(bb):
    num_double_bonds = sum(
        1 for bond in bb.mol.GetBonds()
        if bond.GetBondTypeAsDouble() == 2
    )
    return 100 * num_double_bonds / bb.mol.GetNumBonds()


def _get_functional_group_distance(bb):
    bb = stk.BuildingBlock.init_from_rdkit_mol(
        mol=bb.mol,
        functional_groups=[bb.functional_group],
    )
    return np.mean(list(bb.get_bonder_distances()))


def _get_property_progress(progress, property_fn, bb_filter):
    return [
        np.mean([
            property_fn(bb) for bb in filter(bb_filter, generation)
        ])
        for generation in progress
    ]


def _plot_property(tritopic, ditopic, path):
    label_size = 24
    mpl.rcParams['xtick.labelsize'] = label_size
    mpl.rcParams['ytick.labelsize'] = label_size
    fig, ax = plt.subplots(1, 1, figsize=(23, 5))
    mpl.rc('legend', fontsize=16)

    x = list(range(len(tritopic)))
    ax.plot(
        x,
        ditopic,
        linestyle='--',
        marker='o',
        ms=10,
        color='mediumblue',
        alpha=0.5,
        label='Di-topic',
    )
    ax.plot(
        x,
        tritopic,
        linestyle='--',
        marker='o',
        ms=10,
        color='darkorange',
        alpha=0.5,
        label='Tri-topic',
    )
    ax.legend()
    min_y = max(0, min(val for val in it.chain(tritopic, ditopic))-3)
    max_y = min(100, max(val for val in it.chain(tritopic, ditopic))+3)
    ax.set_ylim([min_y, max_y])

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(path, dpi=1000)
    plt.close()


def _plot_rotatable_bonds(progress, output_directory):
    tritopic = _get_property_progress(
        progress=progress,
        property_fn=_get_percent_rotatable_bonds,
        bb_filter=lambda bb: bb.count == 4,
    )
    ditopic = _get_property_progress(
        progress=progress,
        property_fn=_get_percent_rotatable_bonds,
        bb_filter=lambda bb: bb.count == 6,
    )

    _plot_property(
        tritopic=tritopic,
        ditopic=ditopic,
        path=join(output_directory, 'rotatable_bonds.png'),
    )


def _plot_double_bonds(progress, output_directory):
    tritopic = _get_property_progress(
        progress=progress,
        property_fn=_get_percent_double_bonds,
        bb_filter=lambda bb: bb.count == 4,
    )
    ditopic = _get_property_progress(
        progress=progress,
        property_fn=_get_percent_double_bonds,
        bb_filter=lambda bb: bb.count == 6,
    )

    _plot_property(
        tritopic=tritopic,
        ditopic=ditopic,
        path=join(output_directory, 'double_bonds'),
    )


def _plot_functional_group_distance(progress, output_directory):
    tritopic = _get_property_progress(
        progress=progress,
        property_fn=_get_functional_group_distance,
        bb_filter=lambda bb: bb.count == 4,
    )
    ditopic = _get_property_progress(
        progress=progress,
        property_fn=_get_functional_group_distance,
        bb_filter=lambda bb: bb.count == 6,
    )

    _plot_property(
        tritopic=tritopic,
        ditopic=ditopic,
        path=join(output_directory, 'functional_group_distance.png'),
    )


def main():
    parser = argparse.ArgumentParser()
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

    # Get all the building blocks in a generation.
    progress = [
        [bb for cage in generation for bb in cage.building_blocks]
        for generation in progress
    ]

    for generation in progress:
        for building_block in generation:
            rdkit.SanitizeMol(building_block.mol)
            rdkit.Kekulize(building_block.mol)

    _plot_rotatable_bonds(progress, args.output_directory)
    _plot_double_bonds(progress, args.output_directory)
    _plot_functional_group_distance(progress, args.output_directory)


if __name__ == '__main__':
    main()
