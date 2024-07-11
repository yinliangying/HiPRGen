import os

from HiPRGen.species_filter import species_filter
from HiPRGen.species_questions import (
    li_species_decision_tree,
    positive_penalty,
    species_default_true, no_species_decision_tree
)
from monty.serialization import loadfn, dumpfn
import pickle
import json
from HiPRGen.report_generator import visualize_molecules
from HiPRGen.reaction_filter_payloads import (
    DispatcherPayload,
    WorkerPayload
)
from HiPRGen.bucketing import bucket
from HiPRGen.constants import ROOM_TEMP, Terminal
import subprocess
from HiPRGen.reaction_questions import default_reaction_decision_tree

# def apply_species_filter(pickle_path: str):
#     database_entries = loadfn(pickle_path)
#     # print(database_entries[0])
#     # a = database_entries[0]
#     # print(a.molecule.properties)
#
#     mol_entries = species_filter(
#         database_entries,
#         mol_entries_pickle_location='mol_entries.pickle',
#         species_report='unfiltered_species_report.tex',
#         species_decision_tree=li_species_decision_tree,
#         coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy),
#         generate_unfiltered_mol_pictures=False
#     )
#
# apply_species_filter(r'dump.pickle')


def apply_species_filter(json_path: str):
    with open(json_path, 'r') as f:
        database_entries = json.load(fp=f)

    mol_entries = species_filter(
        database_entries,
        mol_entries_pickle_location='mol_entries.pickle',
        species_report='unfiltered_species_report.tex',
        species_decision_tree=no_species_decision_tree,
        coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy),
        generate_unfiltered_mol_pictures=True
    )


def get_rn_db(pickle_path: str):
    for a_file in os.listdir('./'):
        if a_file.endswith('sqlite'):
            os.remove(a_file)

    with open(pickle_path, 'rb') as file:
        mol_entries = pickle.load(file)
    # mol_entries = loadfn(pickle_path)
    bucket(mol_entries, 'buckets.sqlite')
    params = {
        'temperature': ROOM_TEMP,
        'electron_free_energy': -1.4
    }

    dispatcher_payload = DispatcherPayload(
        'buckets.sqlite',
        'rn.sqlite',
        'reaction_report.tex'
    )

    worker_payload = WorkerPayload(
        'buckets.sqlite',
        default_reaction_decision_tree,
        params,
        Terminal.DISCARD
    )

    dumpfn(dispatcher_payload, 'dispatcher_payload.json')
    dumpfn(worker_payload, 'worker_payload.json')
    subprocess.run(
        [
            'mpirun',
            '--use-hwthread-cpus',
            '-n',
            str(64),
            'python',
            'run_network_generation.py',
            'mol_entries.pickle',
            'dispatcher_payload.json',
            'worker_payload.json'
        ]
    )


apply_species_filter(r'real_dump.json')
get_rn_db('mol_entries.pickle')




