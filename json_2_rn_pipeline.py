import os
import sqlite3
import shutil
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




def apply_species_filter(json_path: str,network_folder: str):
    with open(json_path, 'r') as f:
        database_entries = json.load(fp=f)

    mol_entries = species_filter(
        database_entries,
        mol_entries_pickle_location=f'{network_folder}/mol_entries.pickle',
        species_report=f'{network_folder}/unfiltered_species_report.tex',
        species_decision_tree=no_species_decision_tree,
        coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy),
        generate_unfiltered_mol_pictures=False
    )
    return mol_entries


def get_rn_db(pickle_path: str,network_folder: str):
    # for a_file in os.listdir('../../../../Downloads/'):
    #     if a_file.endswith('sqlite'):
    #         os.remove(a_file)

    with open(pickle_path, 'rb') as file:
        mol_entries = pickle.load(file)
    # mol_entries = loadfn(pickle_path)
    bucket(mol_entries, f'{network_folder}/buckets.sqlite')
    params = {
        'temperature': ROOM_TEMP,
        'electron_free_energy': -1.4
    }

    dispatcher_payload = DispatcherPayload(
        f'{network_folder}/buckets.sqlite',
        f'{network_folder}/rn.sqlite',
        f'{network_folder}/reaction_report.tex'
    )

    worker_payload = WorkerPayload(
        f'{network_folder}/buckets.sqlite',
        default_reaction_decision_tree,
        params,
        Terminal.DISCARD
    )

    dumpfn(dispatcher_payload, f'{network_folder}/dispatcher_payload.json')
    dumpfn(worker_payload, f'{network_folder}/worker_payload.json')
    subprocess.run(
        [
            'mpirun',
            '--use-hwthread-cpus',
            '-n',
            str(64),
            'python',
            'run_network_generation.py',
            pickle_path,
            f'{network_folder}/dispatcher_payload.json',
            f'{network_folder}/worker_payload.json'
        ]
    )

def unit_lib(united_lib_path, new_lib_json_path,old_lib_path):

    #合并old new mol_entries
    new_mol_entries=apply_species_filter(new_lib_path, new_lib_json_path)
    with open(f"{old_lib_path}/mol_entries.pickle", 'rb') as file:
        old_mol_entries = pickle.load(file)

    united_mol_entries=[]
    for mol in old_mol_entries:
        mol.source="old"
        united_mol_entries.append(mol)
    for mol in new_mol_entries:
        mol.source="new"
        mol.ind=len(united_mol_entries)
        united_mol_entries.append(mol)

    with open(f"{united_lib_path}/mol_entries.pickle", 'wb') as f:
        pickle.dump(united_mol_entries, f)

    #在old 基础上建立united rn
    shutil.copy(f"{old_lib_path}/rn.sqlite", f"{united_lib_path}/rn.sqlite")
    rn_con = sqlite3.connect(  f"{united_lib_path}/rn.sqlite")
    rn_cur = rn_con.cursor()
    rn_cur.execute("DROP TABLE IF EXISTS metadata;")
    rn_con.commit()

    bucket(united_mol_entries, f'{united_lib_path}/buckets.sqlite')
    params = {
        'temperature': ROOM_TEMP,
        'electron_free_energy': -1.4
    }
    dispatcher_payload = DispatcherPayload(
        f'{united_lib_path}/buckets.sqlite',
        f'{united_lib_path}/rn.sqlite',
        f'{united_lib_path}/reaction_report.tex'
    )
    worker_payload = WorkerPayload(
        f'{united_lib_path}/buckets.sqlite',
        default_reaction_decision_tree,
        params,
        Terminal.DISCARD
    )
    dumpfn(dispatcher_payload, f'{united_lib_path}/dispatcher_payload.json')
    dumpfn(worker_payload, f'{united_lib_path}/worker_payload.json')
    subprocess.run(
        [
            'mpirun',
            '--use-hwthread-cpus',
            '-n',
            str(64),
            'python',
            'run_network_generation_united.py',
            f"{united_lib_path}/mol_entries.pickle" 
            f'{united_lib_path}/dispatcher_payload.json' 
            f'{united_lib_path}/worker_payload.json'
        ]
    )

if __name__ == '__main__':

    old_lib_path="old_lib"  # old_lib/rn.sqlite  old_lib/buckets.sqlite old_lib/mol_entries.pickle
    new_lib_path="new_lib"
    united_lib_path="united_lib"

    new_lib_json_path=f"{new_lib_path}/real_dump.json"
    # apply_species_filter(new_lib_path, new_lib_json_path)
    # get_rn_db(new_lib_json_path ,new_lib_path)


    unit_lib(united_lib_path, new_lib_json_path,old_lib_path)



