import os
import shutil
import sqlite3
from rdkit import Chem
from rdkit.Chem import AllChem
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge
import pickle
from HiPRGen.report_generator import ReportGenerator
import json
import sys
from HiPRGen.get_common_sub_mol_info import get_common_sub_mol_info
from HiPRGen.launching_entry import libe_default_paths
sys.path.append(r'/usr/local/lib/python3.10/dist-packages')

from pdf2image import convert_from_path

#获取常见辅料分子id  H+ OH- Li2CO3这些

def smiles_to_xyz(smiles: str, output_file_path: str):
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    # Add hydrogen atoms
    mol = Chem.AddHs(mol)

    # Generate a 3D conformation
    AllChem.EmbedMolecule(mol)

    # Optimize the geometry
    AllChem.MMFFOptimizeMolecule(mol)

    # Get the atomic positions
    conf = mol.GetConformer()
    atoms = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]

    # Write the XYZ file
    xyz_content = f'{mol.GetNumAtoms()}\n\n'
    for atom, pos in zip(mol.GetAtoms(), atoms):
        symbol = atom.GetSymbol()
        x, y, z = pos
        xyz_content += f'{symbol} {x:10.6f} {y:10.6f} {z:10.6f}\n'
    with open(output_file_path, 'w') as f:
        f.write(xyz_content)


def tex_to_results(tex_file_path):
    cwd_tex = os.getcwd()
    abs_tex_file_path = os.path.abspath(tex_file_path)
    workbase, tex_file_name = os.path.split(abs_tex_file_path)
    os.chdir(workbase)
    os.system(f'pdflatex {tex_file_name}')
    for a_file in os.listdir('./'):
        if a_file.endswith('.tex') or a_file.endswith('.aux') or a_file.endswith('.log'):
            os.remove(a_file)
    pure_file_name = os.path.splitext(tex_file_name)[0]
    images = convert_from_path(f'{pure_file_name}.pdf')
    png_folder_name = pure_file_name + '_png'
    os.makedirs(png_folder_name)
    os.chdir(png_folder_name)
    for i, image in enumerate(images):
        image.save(f"page_{i + 1}.png", "PNG")
    os.chdir(cwd_tex)



def dump_row_info_as_txt(row_info, txt_file_path, sub_id_file_path, sub_id_list: list):
    with open(txt_file_path, 'w') as f:
        f.write('reaction_id | number_of_reactants | number_of_products | reactant_1 | reactant_2 | product_1 | product_2 | rate | dG | dG_barrier | is_redox\n')
        f.write('---------------------------------------------------------------------------------------------------------------------------------------------\n')
        for row in row_info:
            f.write(f"{row[0]} | {row[1]} | {row[2]} | {row[3]} | {row[4]} | {row[5]} | {row[6]} | {row[7]} | {row[8]} | {row[9]} | {row[10]}")
            f.write('\n')
    new_rev_dict = {}
    for an_id in sub_id_list:
        new_rev_dict.update({an_id: rev_common_sub_mol_info[an_id]})
    with open(sub_id_file_path, 'w') as json_file:
        json.dump(new_rev_dict, json_file)

def main_query(reaction_db_file_path: str, mol_entry_file_path: str, main_mol_smiles: str, mol_picture_folder_path: str):
    with open(mol_entry_file_path, "rb") as pickle_file:
        mol_entries = pickle.load(pickle_file)
    reaction_db_file_path = os.path.abspath(reaction_db_file_path)
    mol_picture_folder_abs_path = os.path.abspath(mol_picture_folder_path)

    cwd_ = os.getcwd()
    if os.path.exists('query_results'):
        shutil.rmtree("query_results")
    os.makedirs('query_results')
    os.chdir('query_results')
    try:
        smiles_to_xyz(smiles=main_mol_smiles, output_file_path=r'main_mol.xyz')
        q_id = find_mol_entry_from_xyz_and_charge(mol_entries=mol_entries,
                                                  xyz_file_path=r'main_mol.xyz',
                                                  charge=0)
    except:
        print('The main molecule has not been found in species database.')
        with open('No species match with the input SMILES', 'w') as f:
            pass
        os.chdir(cwd_)
        return
    db_conn = sqlite3.connect(reaction_db_file_path)
    # Main molecule as reactant
    cursor = db_conn.cursor()
    cursor.execute("SELECT * FROM reactions WHERE reactant_1 = ? OR reactant_2 = ?", (q_id, q_id))
    rows = cursor.fetchall()
    rows.sort(key=lambda row: row[8])

    if not rows:
        print(f"No reactions found with reactant_id: {q_id}, which is the main molecule queried.")
    else:
        main_mol_as_reactant_rows = []
        sub_id_list = []
        for row in rows:
            if row[3] == q_id:
                sub_id = row[4]
            else:
                sub_id = row[3]
            if sub_id in common_sub_mol_ids:
                main_mol_as_reactant_rows.append(row)
                sub_id_list.append(sub_id)
        if len(main_mol_as_reactant_rows) > 0:
            dump_row_info_as_txt(row_info=main_mol_as_reactant_rows, txt_file_path=r'main_mol_as_reactant_info.txt', sub_id_list=sub_id_list, sub_id_file_path=r'reactant_sub_id.json')
            reporter = ReportGenerator(mol_entries=mol_entries,
                                       report_file_path=r'main_mol_as_reactant_info.tex',
                                       fixed_mol_pictures_folder=mol_picture_folder_abs_path)
            reporter.emit_text("Main molecule as reactant report")
            for a_row in main_mol_as_reactant_rows:
                reaction_info = dict(zip(reaction_head, a_row))
                reaction_info.update({'reactants': a_row[3:5], 'products': a_row[5:7]})
                reporter.emit_reaction(reaction_info)
                reporter.emit_newline()
            reporter.finished()
            tex_to_results(tex_file_path='main_mol_as_reactant_info.tex')

    # Main molecule as product
    cursor.execute("SELECT * FROM reactions WHERE product_1 = ? OR product_2 = ?", (q_id, q_id))
    rows = cursor.fetchall()
    rows.sort(key=lambda row: row[8])

    if not rows:
        print(f"No reactions found with product_id: {q_id}, which is the main molecule queried.")
    else:
        main_mol_as_product_rows = []
        sub_id_list = []
        for row in rows:
            if row[5] == q_id:
                sub_id = row[6]
            else:
                sub_id = row[5]
            if sub_id in common_sub_mol_ids:
                main_mol_as_product_rows.append(row)
                sub_id_list.append(sub_id)
        if len(main_mol_as_product_rows) > 0:
            dump_row_info_as_txt(row_info=main_mol_as_product_rows, txt_file_path=r'main_mol_as_product_info.txt', sub_id_list=sub_id_list, sub_id_file_path=r'product_sub_id.json')
            reporter = ReportGenerator(mol_entries=mol_entries,
                                       report_file_path=r'main_mol_as_product_info.tex',
                                       fixed_mol_pictures_folder=mol_picture_folder_abs_path)
            reporter.emit_text("Main molecule as product report")
            for a_row in main_mol_as_product_rows:
                reaction_info = dict(zip(reaction_head, a_row))
                reaction_info.update({'reactants': a_row[3:5], 'products': a_row[5:7]})
                reporter.emit_reaction(reaction_info)
                reporter.emit_newline()
            reporter.finished()
            tex_to_results(tex_file_path='main_mol_as_product_info.tex')
    main_mol_id_filename = f'main_mol_id_is_{q_id}'
    with open(main_mol_id_filename,  'w') as f:
        pass
    os.chdir(cwd_)
    db_conn.close()



with open(r"./app_param.json") as f:
    data = json.load(f)

if "database_dir" not in data:
    database_paths=libe_default_paths
else:
    database_dir = data["database_dir"]
    database_paths = {
        'rn_db_path':os.path.join(database_dir, r'rn.sqlite'),
        'mol_entry_file_path':os.path.join(database_dir, r'mol_entries.pickle'),
        'mol_picture_folder_path':os.path.join(database_dir, r'mol_pictures'),
    }
##################################################################################################################

##################################################################################################################
common_sub_mol_info = get_common_sub_mol_info(
    mol_pkl_path=database_paths['mol_entry_file_path'],
)
rev_common_sub_mol_info = dict(zip(common_sub_mol_info.values(), common_sub_mol_info.keys()))
common_sub_mol_ids = list(common_sub_mol_info.values())
reaction_head = ['reaction_id', 'number_of_reactants', 'number_of_products', 'reactant_1', 'reactant_2', 'product_1', 'product_2', 'rate', 'dG', 'dG_barrier', 'is_redox']

# main_query(
#     reaction_db_file_path=r'/home/postgres/workbase/dpdispatcher/hiprgen_workbase/hiprgen_data/libe/rn.sqlite',
#     mol_entry_file_path=r'/home/postgres/workbase/dpdispatcher/hiprgen_workbase/hiprgen_data/libe/mol_entries.pickle',
#     mol_picture_folder_path=r'/home/postgres/workbase/dpdispatcher/hiprgen_workbase/hiprgen_data/mol_pictures',
#     main_mol_smiles=data['smiles'],
# )


main_query(
    reaction_db_file_path=database_paths['rn_db_path'],
    mol_entry_file_path=database_paths['mol_entry_file_path'],
    mol_picture_folder_path=database_paths['mol_picture_folder_path'],
    main_mol_smiles=data['smiles'],
)
