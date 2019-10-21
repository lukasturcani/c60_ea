import json
import rdkit.Chem.AllChem as rdkit


class BuildingBlockData:
    __slots__ = ['mol', 'functional_group', 'count']

    def __init__(self, mol, functional_group, count):
        self.mol = mol
        self.functional_group = functional_group
        self.count = count


class CageData:
    __slots__ = ['building_blocks', 'cage', 'name', 'topology']

    def __init__(self, building_blocks, cage, name, topology):
        self.building_blocks = building_blocks
        self.cage = cage
        self.name = name
        self.topology = topology

    def __repr__(self):
        return f'Cage({self.name})'

    def __str__(self):
        return repr(self)


def _get_building_blocks(cage_data):
    return tuple(
        BuildingBlockData(
            mol=_get_mol(bb_data['mol_block']),
            functional_group=bb_data['func_grp'],
            count=count,
        )
        for bb_data, count in cage_data['bb_counter']
    )


def _get_mol(mol_block):
    return rdkit.MolFromMolBlock(
        molBlock=mol_block,
        sanitize=False,
        removeHs=False,
    )


def get_database(database_path):
    with open(database_path, 'r') as f:
        database = json.load(f)

    for cage_data in database:
        yield CageData(
            building_blocks=_get_building_blocks(cage_data),
            cage=_get_mol(cage_data['mol_block']),
            name=cage_data['name'],
            topology=cage_data['topology'].split('(')[0],
        )
