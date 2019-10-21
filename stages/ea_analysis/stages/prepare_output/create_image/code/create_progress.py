from utilities import get_database
import argparse
from os.path import join
import pickle


def _get_generations(log):
    names = []
    for line in log:
        if line.isspace():
            yield [name for fitness_value, name in sorted(names)]
            names = []
        else:
            name, *_, fitness_value = line.split()
            names.append((fitness_value, name))


def get_progress(log_path, database_path):
    cages = {cage.name: cage for cage in get_database(database_path)}

    with open(log_path, 'r') as f:
        log = f.readlines()

    for names in _get_generations(log):
        yield [cages[name] for name in names]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'output_folder',
        help=(
            'The folder holding the output of the EA run. It will '
            'hold the files "progress.log" and "database.json" '
            'directly.'
        ),
    )
    parser.add_argument(
        'dump_path',
        help=(
            'The path to the file into which the pickled progress '
            'is dumped.'
        ),
    )

    args = parser.parse_args()
    log_path = join(args.output_folder, 'progress.log')
    database_path = join(args.output_folder, 'database.json')

    progress = list(get_progress(log_path, database_path))

    with open(args.dump_path, 'wb') as f:
        pickle.dump(progress, f)


if __name__ == '__main__':
    main()
