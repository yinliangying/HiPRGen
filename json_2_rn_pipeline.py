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
    print("删除旧文件")
    if os.path.exists(f'{network_folder}/rn.sqlite'):
        os.remove(f'{network_folder}/rn.sqlite')
    if os.path.exists(f'{network_folder}/buckets.sqlite'):
        os.remove(f'{network_folder}/buckets.sqlite')

    # 编译cpp
    os.system("rm /root/HiPRGen/HiPRGen/fragment_matching_found.so")
    print("rm /root/HiPRGen/HiPRGen/fragment_matching_found.so OK")
    os.system("g++ -shared  -O3  -fPIC /root/HiPRGen/HiPRGen/fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so")
    print("g++ -shared  -O3  -fPIC /root/HiPRGen/HiPRGen/fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so OK")

    with open(pickle_path, 'rb') as file:
        mol_entries = pickle.load(file)


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
    number_of_threads = os.popen("nproc").read().strip()
    print(f"number_of_threads:{number_of_threads}")
    subprocess.run(
        [
            'mpirun',
            '--use-hwthread-cpus',
            '-n',
            number_of_threads,
            'python',
            'run_network_generation.py',
            pickle_path,
            f'{network_folder}/dispatcher_payload.json',
            f'{network_folder}/worker_payload.json'
        ]
    )

def unit_rn_db(united_network_folder, new_lib_json_path, new_network_folder, old_network_folder):
    print("删除旧文件")
    if os.path.exists(f"{united_network_folder}/buckets.sqlite"):
        os.remove(f"{united_network_folder}/buckets.sqlite")
    if os.path.exists(f"{united_network_folder}/rn.sqlite"):
        os.remove(f"{united_network_folder}/rn.sqlite")
    if os.path.exists(f"{united_network_folder}/mol_entries.pickle"):
        os.remove(f"{united_network_folder}/mol_entries.pickle")

    print("合并old new mol_entries")
    new_mol_entries=apply_species_filter(new_lib_json_path, new_network_folder)
    with open(f"{old_network_folder}/mol_entries.pickle", 'rb') as file:
        old_mol_entries = pickle.load(file)
    united_mol_entries=[]
    for mol in old_mol_entries:
        mol.source="old"
        united_mol_entries.append(mol)
    for mol in new_mol_entries:
        mol.source="new"
        mol.ind=len(united_mol_entries)
        united_mol_entries.append(mol)
    if not os.path.exists(united_network_folder):
        os.mkdir(united_network_folder)

    #编译cpp
    os.system("rm /root/HiPRGen/HiPRGen/fragment_matching_found.so")
    print("rm /root/HiPRGen/HiPRGen/fragment_matching_found.so OK")
    os.system("g++ -shared  -O3  -fPIC /root/HiPRGen/HiPRGen/fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so")
    print("g++ -shared  -O3  -fPIC /root/HiPRGen/HiPRGen/fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so OK")

    print("保存mol_entries")
    with open(f"{united_network_folder}/mol_entries.pickle", 'wb') as f:
        pickle.dump(united_mol_entries, f)

    print("在old 基础上建立united rn")
    shutil.copy(f"{old_network_folder}/rn.sqlite", f"{united_network_folder}/rn.sqlite")
    rn_con = sqlite3.connect(  f"{united_network_folder}/rn.sqlite")
    rn_cur = rn_con.cursor()
    rn_cur.execute("DROP TABLE IF EXISTS metadata;")
    rn_con.commit()


    bucket(united_mol_entries, f'{united_network_folder}/buckets.sqlite')
    params = {
        'temperature': ROOM_TEMP,
        'electron_free_energy': -1.4
    }
    dispatcher_payload = DispatcherPayload(
        f'{united_network_folder}/buckets.sqlite',
        f'{united_network_folder}/rn.sqlite',
        f'{united_network_folder}/reaction_report.tex'
    )
    worker_payload = WorkerPayload(
        f'{united_network_folder}/buckets.sqlite',
        default_reaction_decision_tree,
        params,
        Terminal.DISCARD
    )
    dumpfn(dispatcher_payload, f'{united_network_folder}/dispatcher_payload.json')
    dumpfn(worker_payload, f'{united_network_folder}/worker_payload.json')
    number_of_threads = os.popen("nproc").read().strip()
    print(f"number_of_threads:{number_of_threads}")
    subprocess.run(
        [
            'mpirun',
            '--use-hwthread-cpus',
            '-n',
            "2",#number_of_threads,
            'python',
            'run_network_generation.py',
            f"{united_network_folder}/mol_entries.pickle",
            f'{united_network_folder}/dispatcher_payload.json',
            f'{united_network_folder}/worker_payload.json'
        ]
    )

if __name__ == '__main__':
    # # workdir /root/test_fmol_unfilter    ~/HiPRGen/HiPRGen# cp  fragment_matching_found.cpp  fragment_matching_found_cpp.py 以及这个文件   /root/test_fmol_unfilter/

    old_network_folder= "old_lib"  # old_lib/rn.sqlite  old_lib/buckets.sqlite old_lib/mol_entries.pickle
    new_network_folder= "new_lib"
    united_network_folder= "united_lib"
    new_lib_json_path=f"{new_network_folder}/real_dump.json"
    apply_species_filter(new_lib_json_path, new_network_folder)
    get_rn_db("new_lib/mol_entries.pickle", new_network_folder)


    unit_rn_db(united_network_folder, new_lib_json_path, new_network_folder, old_network_folder)



