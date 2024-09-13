import sqlite3
reaction_db_file_path=rf'/root/HiPRGen/data/new_libe_fmol_20240731/rn.sqlite'
db_conn = sqlite3.connect(reaction_db_file_path)
cursor = db_conn.cursor()
q_id="142"
cursor.execute("SELECT * FROM reactions WHERE reactant_1 = ? OR reactant_2 = ?", (q_id, q_id))
rows = cursor.fetchall()
for row in rows:
    print(row)
import HiPRGen

from rdkit import Chem
from rdkit.Chem import AllChem
import pickle
mol_pickle_file="/root/HiPRGen/data/new_libe_fmol_20240731/mol_entries.pickle"
mol_pickle=pickle.load(open(mol_pickle_file,'rb'))
for mol_entry in mol_pickle:
    molecule = mol_entry.molecule
    molecule.to(filename="mol_path")
    mol = AllChem.MolFromMolFile("mol_path")
    try:
        Chem.SanitizeMol(mol)
        AllChem.RemoveStereochemistry(mol)
        smiles = AllChem.MolToSmiles(mol)
    except:
        smiles = ""
    print(smiles)