import os
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge
import pickle
import sys

sys.path.append(r'/usr/local/lib/python3.10/dist-packages')
sub_mol_xyz_path = r'/root/test_common_xyz/common_xyz_files'
def get_common_sub_mol_info(mol_pkl_path: str,):
    cwd_sub = os.getcwd()
    with open(mol_pkl_path, "rb") as pickle_file:
        mol_entries = pickle.load(pickle_file)
    q_id_list = []
    mol_name_list = []
    for a_charge_folder in os.listdir(sub_mol_xyz_path):
        os.chdir(sub_mol_xyz_path)
        os.chdir(a_charge_folder)
        if a_charge_folder.startswith('zero'):
            a_charge = 0
        elif a_charge_folder.startswith('positive_one_'):
            a_charge = 1
        else:
            a_charge = -1

        for a_file in os.listdir('./'):
            try:
                q_id = find_mol_entry_from_xyz_and_charge(mol_entries=mol_entries,
                                                          xyz_file_path=a_file,
                                                          charge=a_charge)
                q_id_list.append(q_id)
                mol_name_list.append(a_file[:-4])
            except:
                pass
    common_sub_mol_info = dict(zip(mol_name_list, q_id_list))
    os.chdir(cwd_sub)
    return common_sub_mol_info
