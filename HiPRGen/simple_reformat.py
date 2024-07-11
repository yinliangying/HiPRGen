from ase.db import connect
from ase.atoms import Atoms
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.core.structure import Molecule
import json
import numpy as np
from monty.serialization import loadfn, dumpfn


def convert_ndarray_to_list(obj):
    if isinstance(obj, dict):
        return {k: convert_ndarray_to_list(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_ndarray_to_list(elem) for elem in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

def transform_db_2_entry(db_path: str, pickle_path: str):
    mol_graphs = []
    with connect(db_path) as db:
        for idx, row in enumerate(db.select()):
            data = row.data
            data['partial_charges'] = {
                'mulliken': data.atomcharges['mulliken'].tolist(),
                'nbo': data.atomcharges['natural'].tolist(),
            }
            del data['atomcharges']
            # data = convert_ndarray_to_list(data)
            atoms = row.toatoms()
            symbols = [atom.symbol for atom in atoms]
            positions = [atom.position.tolist() for atom in atoms]
            molecule = Molecule(symbols, positions, charge=data['charge'], spin_multiplicity=data['mult'], properties=data)
            mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())
            final_info = mol_graph.as_dict()
            final_info['graphs']['nodes'] = convert_ndarray_to_list(final_info['graphs']['nodes'])
            final_info['xyz'] = positions
            final_info['partial_spins'] = {
                'mulliken': None,
                'nbo': None
            }
            # final_info = convert_ndarray_to_list(final_info)
            mol_graphs.append(final_info)

            if len(mol_graphs) == 50:
                break

    dumpfn(mol_graphs, pickle_path)

def transform_db_2_json(db_path: str, json_path: str):
    mol_graphs = []
    with connect(db_path) as db:
        for idx, row in enumerate(db.select()):
            data = row.data
            data['partial_charges'] = {
                'mulliken': data.atomcharges['mulliken'].tolist(),
                'nbo': data.atomcharges['natural'].tolist(),
            }
            del data['atomcharges']
            # data = convert_ndarray_to_list(data)
            atoms = row.toatoms()
            symbols = [atom.symbol for atom in atoms]
            positions = [atom.position.tolist() for atom in atoms]
            molecule = Molecule(symbols, positions, charge=data['charge'], spin_multiplicity=data['mult'], properties=data)
            mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())
            final_info = mol_graph.as_dict()
            final_info['graphs']['nodes'] = convert_ndarray_to_list(final_info['graphs']['nodes'])
            final_info['xyz'] = positions
            final_info['partial_spins'] = {
                'mulliken': None,
                'nbo': None
            }
            # final_info = convert_ndarray_to_list(final_info)
            mol_graphs.append(final_info)
            # if len(mol_graphs) == 50:
            #     break
    with open(json_path, 'w') as f:
        json.dump(mol_graphs, f)


transform_db_2_json(db_path=r'dump.db', json_path=r'real_dump.json')

# transform_db_2_entry(db_path=r'dump.db', pickle_path=r'dump.pickle')