import os
import shutil
from pathlib import Path

from HiPRGen.launching_entry import run_with_id, run_with_smiles

#
# class Subsidiary_Molecule_Options(String, Enum):
#     """
#     Define the target to be predicted.
#     """
#     Li = 'Li+1'
#     OH = 'OH-1'
#     C2H4 = 'C2H4'
#     H2O = 'H2O'
#     H2 = 'H2'
#     CO = 'CO'
#     CO2 = 'CO2'
#     LiF = 'LiF'
#     nll = 'None'
#     LiCO3 = 'LiCO3-1'
#     H = 'H+1'
#     Li2CO3 = 'Li2CO3+1'



if __name__ == '__main__':
    import sys
    import json
    #python /root/HiPRGen/hiprgen_kmc_bohrium_entry.py kMC_pathfinding --json-config app_param.json
    if len(sys.argv) !=4:
        print('Please run with kMC_pathfinding --json-config app_param.json')
        exit()
    if sys.argv[1] != 'kMC_pathfinding':
        print('Please run with kMC_pathfinding')
        exit()

    if sys.argv[2] != '--json-config':
        print('Please run with --json-config')
        exit()

    json_param_path = sys.argv[3]

    # {
    #     "output_dir": output_dir,
    #     "Input_Format": {
    #         "type": "ID",
    #         "main_mol_id": 13589,
    #         "sub_mol_ids": [3274, 14191, 145],
    #         "n_sim": 1000,
    #         "database_dir": "/root/HiPRGen/data/libe_and_fmol_0911_all_rn_filter",
    #         # "/root/HiPRGen/data/new_libe_fmol_20240731",#"/root/HiPRGen/data/libe",#"/root/HiPRGen/data/new_libe_fmol_20240731",
    #     },
    # }
    param_dict= json.load(open(json_param_path))

    database_dir=param_dict["Input_Format"]["database_dir"]
    main_mol_id=param_dict["Input_Format"]["main_mol_id"]
    sub_mol_ids=param_dict["Input_Format"]["sub_mol_ids"]
    n_sim=param_dict["Input_Format"]["n_sim"]
    output_directory=param_dict["output_dir"]
    database_paths = {
        'rn_db_path': rf'{database_dir}/rn.sqlite',
        'mol_entry_file_path': rf'{database_dir}/mol_entries.pickle',
        'mol_picture_folder_path': rf'{database_dir}/mol_pictures',
    }
    run_with_id(
        main_mol_id= main_mol_id ,
        sub_mol_ids=sub_mol_ids,
        simulation_times=n_sim,
        num_cores=8,
        output_dir=output_directory,
        default_file_paths=database_paths
    )
