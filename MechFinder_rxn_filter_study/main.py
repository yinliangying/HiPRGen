import sys
import os
import argparse
import os
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from ase.db import connect
from ase import Atoms
from ase.symbols import symbols2numbers
from openbabel import openbabel
from HiPRGen.mol_pkl_2_mol_pics import xyz_2_db_mol,obmol_to_rdkit_mol
import pandas as pd
import sqlite3
import os
from tqdm import tqdm

def set_radical_electrons(rd_mol, mol_charge): #mol_charge 0 -1 +1


    star_atom_num=0
    for idx, atom in enumerate(rd_mol.GetAtoms()):
        atom_num = atom.GetAtomicNum()
        if atom_num == 3:
            continue
        typical_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
        actual_valence = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        if actual_valence < typical_valence:
            star_atom_num+=1
            continue

    if (mol_charge==0 and star_atom_num ==1 ) or (abs(mol_charge)==star_atom_num): #这些情况下电荷和和自由基能够无歧义分配
        if mol_charge==0 and star_atom_num ==1:
            for idx, atom in enumerate(rd_mol.GetAtoms()):
                atom_num = atom.GetAtomicNum()
                if atom_num == 3:
                    continue
                typical_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
                actual_valence = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
                if actual_valence < typical_valence:
                    rd_mol.GetAtomWithIdx(idx).SetNumRadicalElectrons(int(typical_valence-actual_valence))
                    break
            return rd_mol,True
        else:
            if (mol_charge == 0 and star_atom_num == 0):
                return rd_mol,True
            else:
                charge_sign = -1 if mol_charge < 0 else 1
                for idx, atom in enumerate(rd_mol.GetAtoms()):
                    atom_num = atom.GetAtomicNum()
                    if atom_num == 3:
                        continue
                    typical_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
                    actual_valence = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
                    if actual_valence < typical_valence:
                        rd_mol.GetAtomWithIdx(idx).SetFormalCharge(charge_sign)
                return rd_mol,True

    else:
        return rd_mol,False
def filter_mol_entries(pickle_path: str,output_smi_csv_path: str) -> str:
    """
    过滤分子，区分哪些分子的缺键原子是完全可标注的（目前只考察这类分子） 大约1/2
    """
    # Load pickle file
    with open(pickle_path, 'rb') as f:
        f_data = pickle.load(f)

    with open(output_smi_csv_path,"w") as fp_out:
        print("idx,smiles,well_define",file=fp_out)
        for idx, a_mol_info in enumerate(f_data):
            # Create XYZ file
            a_mol = a_mol_info.molecule
            a_new_atoms = Atoms(symbols=a_mol.labels, positions=a_mol.cart_coords)

            xyz_filename = f"{idx}.xyz"
            a_new_atoms.write(xyz_filename)

            # Convert XYZ to OBMol
            obmol = xyz_2_db_mol(idx)

            # Convert OBMol to RDKit Mol
            rdkit_mol = obmol_to_rdkit_mol(obmol)

            # Set charge and spin multiplicity
            rdkit_mol.SetProp("_Name", f"Molecule {idx}")
            rdkit_mol.SetProp("Charge", str(a_mol.charge))
            rdkit_mol.SetProp("SpinMultiplicity", str(a_mol.spin_multiplicity))
            rdkit_mol,well_define = set_radical_electrons(rdkit_mol, a_mol.charge)

            # Remove temporary XYZ file
            os.remove(xyz_filename)

            smiles= Chem.MolToSmiles(rdkit_mol)
            try:
                smiles=Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
                if smiles=="" or None:
                    continue
            except:
                continue
            print(f"{idx},{smiles},{0 if well_define==False else 1}",file=fp_out)



def trans_rxn_db2smarts(smi_csv_path: str,rn_db_path: str,rxn_smarts_output_file: str):
    """
    rn_db_path  exp:"/root/HiPRGen/data/libe_and_fmol_0911_all/rn.sqlite"
    rxn_smarts_output_file = "/personal/Bohrium_task_hiprgen_rn/hiprgen_json2rn_output/libe_and_fmol_0911_all/rxn_smarts.txt"

    """


    id_smiles_dict = {}
    smiles_id_dict = {}
    smi_df= pd.read_csv(smi_csv_path)
    for i, row in tqdm(smi_df.iterrows(),total=smi_df.shape[0]):
        mol_id = int(row["idx"])
        smiles = row["smiles"]
        id_smiles_dict[mol_id] = smiles
        smiles_id_dict[smiles] = mol_id


    rn_con = sqlite3.connect(rn_db_path)
    rn_cur = rn_con.cursor()

    sql_limit = 1000000
    rn_cur.execute(
        f"select  reaction_id, number_of_reactants, number_of_products, reactant_1, reactant_2, product_1, product_2 from "
        f"reactions where reaction_id%1000=0 limit {sql_limit}")

    with open(rxn_smarts_output_file, "w") as fp_out:
        print("reaction_id,rxn_smarts",file=fp_out)
        for row in tqdm(rn_cur, total=sql_limit):
            reaction_id = row[0]
            number_of_reactants = int(row[1])
            number_of_products = int(row[2])

            reactant_1 = int(row[3])
            reactant_2 = int(row[4])
            product_1 = int(row[5])
            product_2 = int(row[6])

            try:
                reactant_1_smiles = id_smiles_dict[reactant_1]
            except:
                continue
            if number_of_reactants == 2:
                try:
                    reactant_2_smiles = id_smiles_dict[reactant_2]
                except:
                    continue
            else:
                reactant_2_smiles = ""
            try:
                product_1_smiles = id_smiles_dict[product_1]
            except:
                continue
            if number_of_products == 2:
                try:
                    product_2_smiles = id_smiles_dict[product_2]
                except:
                    continue
            else:
                product_2_smiles = ""
            if number_of_reactants == 2 and number_of_products == 2:
                rxn_smarts = f"{reactant_1_smiles}.{reactant_2_smiles}>>{product_1_smiles}.{product_2_smiles}"
            elif number_of_reactants == 2 and number_of_products == 1:
                rxn_smarts = f"{reactant_1_smiles}.{reactant_2_smiles}>>{product_1_smiles}"
            elif number_of_reactants == 1 and number_of_products == 2:
                rxn_smarts = f"{reactant_1_smiles}>>{product_1_smiles}.{product_2_smiles}"
            elif number_of_reactants == 1 and number_of_products == 1:
                rxn_smarts = f"{reactant_1_smiles}>>{product_1_smiles}"
            else:
                print(f"error:{number_of_reactants},{number_of_products} {reaction_id}")
                continue

            print(f"{reaction_id},{rxn_smarts}",file=fp_out)




def mapping_rxn(rxn_smarts_file: str, mapped_rxn_smarts_output_file: str):
    from rxnmapper import RXNMapper
    rxn_mapper = RXNMapper()
    rxn_df = pd.read_csv(rxn_smarts_file)
    batch_size=32
    rxn_smarts_batch_list=[]
    rxn_id_batch_list=[]

    with open( mapped_rxn_smarts_output_file, "w") as fp_out:
        print("reaction_id,mapped_rxn",file=fp_out)
        for i, row in tqdm(rxn_df.iterrows(),total=rxn_df.shape[0]):
            rxn_smarts = row["rxn_smarts"]
            reaction_id= row["reaction_id"]


            rxn_id_batch_list.append(reaction_id)
            rxn_smarts_batch_list.append(rxn_smarts)

            if len(rxn_id_batch_list)==batch_size:
                results = rxn_mapper.get_attention_guided_atom_maps(rxn_smarts_batch_list)
                for i, result in enumerate(results):
                    rxn_id=rxn_id_batch_list[i]
                    rxn_smarts=rxn_smarts_batch_list[i]
                    try:
                        mapped_rxn=result["mapped_rxn"]
                        confidence=result["confidence"]
                    except:
                        print(f"error:{rxn_id}")
                        continue
                    print(f"{rxn_id},{mapped_rxn}",file=fp_out)

                rxn_id_batch_list=[]
                rxn_smarts_batch_list=[]

        results = rxn_mapper.get_attention_guided_atom_maps(rxn_smarts_batch_list)
        for i, result in enumerate(results):
            rxn_id = rxn_id_batch_list[i]
            rxn_smarts = rxn_smarts_batch_list[i]
            try:
                mapped_rxn = result["mapped_rxn"]
                confidence = result["confidence"]
            except:
                print(f"error:{rxn_id}")
                continue
            print(f"{rxn_id},{mapped_rxn}", file=fp_out)



def apply_MechFinder(mapped_rxn_smarts_file: str, mech_output_file: str):
    from MechFinder import MechFinder
    finder = MechFinder(collection_dir='MechFinder/collections')

    df=pd.read_csv(mapped_rxn_smarts_file)

    with open(mech_output_file, "w") as fp_out:
        print("reaction_id,rxn_str,updated_reaction,LRT,MT_class,electron_path",file=fp_out)
        for i,row in tqdm(df.iterrows(),total=df.shape[0]):
            rxn_id=row["reaction_id"]
            rxn_str = row["mapped_rxn"]
            try:
                updated_reaction, LRT, MT_class, electron_path = finder.get_electron_path(rxn_str)
            except:
                continue
            if not isinstance(finder.check_exception(MT_class), str):
                if MT_class!="mechanism not in collection":
                    print(f"{rxn_id},{rxn_str},{updated_reaction},{LRT},{MT_class},{electron_path}",file=fp_out)

if __name__ == "__main__":



    data_dir="/personal/Bohrium_task_hiprgen_rn/hiprgen_json2rn_output/libe_and_fmol_0911_all/"
    mol_entries_file=f"{data_dir}mol_entries.pickle"

    filter_mol_entries(mol_entries_file, f"{data_dir}smiles.csv")


    trans_rxn_db2smarts(f"{data_dir}smiles.csv",
                        rn_db_path="/root/HiPRGen/data/libe_and_fmol_0911_all/rn.sqlite",
                        rxn_smarts_output_file=f"{data_dir}rxn_smarts.csv")
    mapping_rxn(f"{data_dir}rxn_smarts.csv",f"{data_dir}rxn_smarts_mapped.csv")
    apply_MechFinder(f"{data_dir}rxn_smarts_mapped.csv",f"{data_dir}rxn_smarts_mapped_mech.csv")

