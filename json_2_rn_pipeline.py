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




def apply_species_filter(json_path: str, output_network_folder: str):
    with open(json_path, 'r') as f:
        database_entries = json.load(fp=f)
    if os.path.exists(f'{output_network_folder}/mol_entries.pickle'):
        with open(f'{output_network_folder}/mol_entries.pickle', 'rb') as file:
            mol_entries = pickle.load(file)
    else:
        mol_entries = species_filter(
            database_entries,
            mol_entries_pickle_location=f'{output_network_folder}/mol_entries.pickle',
            species_report=f'{output_network_folder}/unfiltered_species_report.tex',
            species_decision_tree=no_species_decision_tree,
            coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy),
            generate_unfiltered_mol_pictures=False
        )
    return mol_entries


def get_rn_db(pickle_path: str, output_network_folder: str,machine_num: int,machine_id:int):
    print("删除旧rn.sqlite文件")
    if os.path.exists(f'{output_network_folder}/rn.sqlite'):
        os.remove(f'{output_network_folder}/rn.sqlite')
    # if os.path.exists(f'{output_network_folder}/buckets.sqlite'):
    #     os.remove(f'{output_network_folder}/buckets.sqlite')

    # 编译cpp
    os.system("rm /root/HiPRGen/HiPRGen/fragment_matching_found.so")
    print("rm /root/HiPRGen/HiPRGen/fragment_matching_found.so OK")
    os.system("g++ -shared  -O3  -fPIC /root/HiPRGen/HiPRGen/fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so")
    print("g++ -shared  -O3  -fPIC /root/HiPRGen/HiPRGen/fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so OK")

    with open(pickle_path, 'rb') as file:
        mol_entries = pickle.load(file)

    if not os.path.exists(f'{output_network_folder}/buckets.sqlite'):
        bucket(mol_entries, f'{output_network_folder}/buckets.sqlite')
    params = {
        'temperature': ROOM_TEMP,
        'electron_free_energy': -1.4
    }
    dispatcher_payload = DispatcherPayload(
        f'{output_network_folder}/buckets.sqlite',
        f'{output_network_folder}/rn.sqlite',
        f'{output_network_folder}/reaction_report.tex',
        machine_num=machine_num,machine_id=machine_id,

    )
    worker_payload = WorkerPayload(
        f'{output_network_folder}/buckets.sqlite',
        default_reaction_decision_tree,
        params,
        Terminal.DISCARD
    )
    dumpfn(dispatcher_payload, f'{output_network_folder}/dispatcher_payload.json')
    dumpfn(worker_payload, f'{output_network_folder}/worker_payload.json')
    number_of_threads = os.popen("nproc").read().strip()
    print(f"number_of_threads:{number_of_threads}")
    subprocess.run(
        [
            'mpirun',
            '--use-hwthread-cpus',
            '-n',
            number_of_threads,
            'python',
            '/root/HiPRGen/run_network_generation.py',
            pickle_path,
            f'{output_network_folder}/dispatcher_payload.json',
            f'{output_network_folder}/worker_payload.json'
        ]
    )

def unit_rn_db(output_network_folder, new_lib_json_path, new_network_folder, old_network_folder):
    print("删除旧文件")
    if os.path.exists(f"{output_network_folder}/buckets.sqlite"):
        os.remove(f"{output_network_folder}/buckets.sqlite")
    if os.path.exists(f"{output_network_folder}/rn.sqlite"):
        os.remove(f"{output_network_folder}/rn.sqlite")
    if os.path.exists(f"{output_network_folder}/mol_entries.pickle"):
        os.remove(f"{output_network_folder}/mol_entries.pickle")

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
    if not os.path.exists(output_network_folder):
        os.mkdir(output_network_folder)

    #编译cpp
    os.system("rm /root/HiPRGen/HiPRGen/fragment_matching_found.so")
    print("rm /root/HiPRGen/HiPRGen/fragment_matching_found.so OK")
    os.system("g++ -shared  -O3  -fPIC /root/HiPRGen/HiPRGen/fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so")
    print("g++ -shared  -O3  -fPIC /root/HiPRGen/HiPRGen/fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so OK")

    print("保存mol_entries")
    with open(f"{output_network_folder}/mol_entries.pickle", 'wb') as f:
        pickle.dump(united_mol_entries, f)

    print("在old 基础上建立united rn")
    shutil.copy(f"{old_network_folder}/rn.sqlite", f"{output_network_folder}/rn.sqlite")
    rn_con = sqlite3.connect(  f"{output_network_folder}/rn.sqlite")
    rn_cur = rn_con.cursor()
    rn_cur.execute("DROP TABLE IF EXISTS metadata;")
    rn_con.commit()


    bucket(united_mol_entries, f'{output_network_folder}/buckets.sqlite')
    params = {
        'temperature': ROOM_TEMP,
        'electron_free_energy': -1.4
    }
    dispatcher_payload = DispatcherPayload(
        f'{output_network_folder}/buckets.sqlite',
        f'{output_network_folder}/rn.sqlite',
        f'{output_network_folder}/reaction_report.tex'
    )
    worker_payload = WorkerPayload(
        f'{output_network_folder}/buckets.sqlite',
        default_reaction_decision_tree,
        params,
        Terminal.DISCARD
    )
    dumpfn(dispatcher_payload, f'{output_network_folder}/dispatcher_payload.json')
    dumpfn(worker_payload, f'{output_network_folder}/worker_payload.json')
    number_of_threads = os.popen("nproc").read().strip()
    print(f"number_of_threads:{number_of_threads}")
    subprocess.run(
        [
            'mpirun',
            '--use-hwthread-cpus',
            '-n',
            number_of_threads,
            'python',
            '/root/HiPRGen/run_network_generation.py',
            f"{output_network_folder}/mol_entries.pickle",
            f'{output_network_folder}/dispatcher_payload.json',
            f'{output_network_folder}/worker_payload.json'
        ]
    )

if __name__ == '__main__':
    # test shell
    # cd /personal/Bohrium_task_hiprgen/hiprgen_json2rn_input
    # python /root/HiPRGen/json_2_rn_pipeline.py  -m ab_initio -j new_lib/real_dump.json -o new_lib
    # python /root/HiPRGen/json_2_rn_pipeline.py   -m append -j new_lib/real_dump.json -o united_lib -n old_lib
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mode", help="Mode must be ab_initio or append", type=str, choices=['ab_initio', 'append'], required=True)
    parser.add_argument("-j", "--json_path", help="Path to json file to be calculated", type=str, required=True)
    parser.add_argument("-o", "--output_network_folder", help="Path to output network folder", type=str, required=True)
    parser.add_argument("-n", "--old_network_folder", help="Only for append mode", type=str, required=False)
    parser.add_argument("-a", "--machine_num", help="", type=str, required=False,default="None")
    parser.add_argument("-i", "--machine_id", help="", type=str, required=False,default="None")

    args = parser.parse_args()
    if args.machine_num=="None":
        args.machine_num=None
    else:
        args.machine_num=int(args.machine_num)

    if args.machine_id=="None":
        args.machine_id=None
    else:
        args.machine_id=int(args.machine_id)
    if args.machine_id==None and args.machine_num!=None:
        print("machine_id must be provided if machine_num is provided")
        exit()
    if args.machine_id!=None and args.machine_num==None:
        print("machine_num must be provided if machine_id is provided")
        exit()
    if args.machine_id>=args.machine_num:
        print("machine_id must be smaller than machine_num")
        exit()
    if args.machine_id<0 or args.machine_num<=0:
        print("machine_id and machine_num must be positive")
        exit()



    print(f"Selected mode: {args.mode}")
    if args.mode == 'ab_initio':
        apply_species_filter(args.json_path, args.output_network_folder)
        get_rn_db(f"{args.output_network_folder}/mol_entries.pickle", args.output_network_folder,
                  args.machine_num,args.machine_id)

    if args.mode == 'append':
        new_network_folder= "new_lib"
        unit_rn_db(args.output_network_folder, args.json_path, new_network_folder, args.old_network_folder)



