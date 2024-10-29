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
from rdkit.Chem.Draw import ReactionToImage
import shutil
from PIL import Image, ImageDraw, ImageFont
from localmapper import localmapper

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
            return rd_mol,True,star_atom_num
        else:
            if (mol_charge == 0 and star_atom_num == 0):
                return rd_mol,True,star_atom_num
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
                return rd_mol,True,star_atom_num

    else:
        return rd_mol,False,star_atom_num
def filter_mol_entries(pickle_path: str,output_smi_csv_path: str) -> str:
    """
    过滤分子，区分哪些分子的缺键原子是完全可标注的（目前只考察这类分子） 大约1/2
    """
    # Load pickle file
    with open(pickle_path, 'rb') as f:
        f_data = pickle.load(f)

    with open(output_smi_csv_path,"w") as fp_out:
        print("idx,smiles,well_define,mol_charge,star_atom_num",file=fp_out)
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
            rdkit_mol,well_define,star_atom_num = set_radical_electrons(rdkit_mol, a_mol.charge)

            # Remove temporary XYZ file
            os.remove(xyz_filename)

            smiles= Chem.MolToSmiles(rdkit_mol)
            try:
                smiles=Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
                if smiles=="" or None:
                    continue
            except:
                continue
            print(f"{idx},{smiles},{0 if well_define==False else 1},{a_mol.charge},{star_atom_num}",file=fp_out)

def filter_rxn( smi_csv_path: str,rxn_db_path: str,filtered_rxn_db_path_path: str):
    """
       目前只考虑无自由基无净电荷分子
         2. 只包含通过模板过滤的反应
           3. 如果反应匹配上模板，模板的反应物产物数必须和反应的反应物产物数相同
    """

    # Load smiles
    id_smiles_dict = {}
    smiles_id_dict = {}
    smi_df= pd.read_csv(smi_csv_path)
    for i, row in tqdm(smi_df.iterrows(),total=smi_df.shape[0]):
        mol_id = int(row["idx"])
        smiles = row["smiles"]
        well_define = int(row["well_define"])
        mol_charge = int(row["mol_charge"])
        star_atom_num=int(row["star_atom_num"])
        if well_define==1 and mol_charge==0 and star_atom_num==0:
            id_smiles_dict[mol_id] = smiles
            smiles_id_dict[smiles] = mol_id
    print(f"smiles_id_dict length:{len(smiles_id_dict)}")
    #prepare filtered_rxn_db
    if os.path.exists(filtered_rxn_db_path_path):
        os.remove(filtered_rxn_db_path_path)
    filtered_rxn_con = sqlite3.connect( filtered_rxn_db_path_path)
    filtered_rxn_cur = filtered_rxn_con.cursor()
    filtered_rxn_cur.execute("""
        CREATE TABLE reactions (
                reaction_id         INTEGER NOT NULL PRIMARY KEY,
                number_of_reactants INTEGER NOT NULL,
                number_of_products  INTEGER NOT NULL,
                reactant_1          INTEGER NOT NULL,
                reactant_2          INTEGER NOT NULL,
                product_1           INTEGER NOT NULL,
                product_2           INTEGER NOT NULL,
                rate                REAL NOT NULL,
                dG                  REAL NOT NULL,
                dG_barrier          REAL NOT NULL,
                is_redox            INTEGER NOT NULL,
                template            STRING,
                mapped_rxn          STRING,
                rxn                 STRING
        );
    """)
    filtered_rxn_con.commit()



    rn_con = sqlite3.connect(rxn_db_path)
    rn_cur = rn_con.cursor()

    # rn_cur.execute(f"select count(*) from reactions")
    # sql_limit = rn_cur.fetchone()[0]

    rn_cur.execute(
        f"select  reaction_id, number_of_reactants, number_of_products, reactant_1, reactant_2, product_1, product_2, "
        f"rate,dG,dG_barrier,is_redox from reactions")
    mapping_batch_size=100
    commit_freq=100
    mapping_times = 0
    tmp_mapping_row_list=[]
    tmp_mapping_rxn_list=[]
    mapper = localmapper("cpu")
    for row in tqdm(rn_cur):

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
        try:
            product_1_smiles = id_smiles_dict[product_1]
        except:
            continue
        if number_of_products == 2:
            try:
                product_2_smiles = id_smiles_dict[product_2]
            except:
                continue
        rxn_smiles=reactant_1_smiles
        if number_of_reactants == 2:
            rxn_smiles += "." + reactant_2_smiles
        rxn_smiles += f">>{product_1_smiles}"
        if number_of_products == 2:
            rxn_smiles += "." + product_2_smiles

        tmp_mapping_rxn_list.append(rxn_smiles)
        tmp_mapping_row_list.append(row)

        if len(tmp_mapping_row_list)==mapping_batch_size:
            print("mapping times:",mapping_times)
            tmp_result_list = mapper.get_atom_map(tmp_mapping_rxn_list, return_dict=True)
            print("mapping finished")
            mapping_times += 1
            for tmp_row, tmp_result, tmp_rxn in zip(tmp_mapping_row_list, tmp_result_list,tmp_mapping_rxn_list):
                if tmp_result["confident"] == True:
                    template = tmp_result["template"]
                    mapped_rxn = tmp_result["mapped_rxn"]
                    rxn_reactant, rxn_product = mapped_rxn.split(">>")
                    rxn_reactant_num = len(rxn_reactant.split("."))
                    rxn_product_num = len(rxn_product.split("."))
                    template_reactant, template_product = template.split(">>")
                    template_reactant_num = len(template_reactant.split("."))
                    template_product_num = len(template_product.split("."))
                    if rxn_reactant_num == template_reactant_num and rxn_product_num == template_product_num:
                        filtered_sql_str = f"""
                        insert into reactions 
                        (reaction_id, number_of_reactants,number_of_products,reactant_1,  reactant_2,  product_1,   product_2,  rate,         dG,          dG_barrier,  is_redox,
                        template,mapped_rxn,rxn)
                         values 
                        ({tmp_row[0]},{tmp_row[1]},       {tmp_row[2]},      {tmp_row[3]},{tmp_row[4]},{tmp_row[5]},{tmp_row[6]},{tmp_row[7]},{tmp_row[8]},{tmp_row[9]},{tmp_row[10]},
                        {template},{mapped_rxn},{tmp_rxn})
                        """
                        try:
                            filtered_rxn_cur.execute(filtered_sql_str)
                        except:
                            print("error")
                            print(filtered_rxn_cur)
                            exit(1)
            tmp_mapping_rxn_list=[]
            tmp_mapping_row_list=[]

            if mapping_times % commit_freq == 0:
                print(f"commit at {mapping_times}")
                filtered_rxn_con.commit()

    tmp_result_list = mapper.get_atom_map(tmp_mapping_rxn_list, return_dict=True)
    for tmp_row, tmp_result, tmp_rxn in zip(tmp_mapping_row_list, tmp_result_list, tmp_mapping_rxn_list):
        if tmp_result["confident"] == True:
            template = tmp_result["template"]
            mapped_rxn = tmp_result["mapped_rxn"]
            rxn_reactant, rxn_product = mapped_rxn.split(">>")
            rxn_reactant_num = len(rxn_reactant.split("."))
            rxn_product_num = len(rxn_product.split("."))
            template_reactant, template_product = template.split(">>")
            template_reactant_num = len(template_reactant.split("."))
            template_product_num = len(template_product.split("."))
            if rxn_reactant_num == template_reactant_num and rxn_product_num == template_product_num:
                filtered_sql_str = f"""
                            insert into reactions 
                            (reaction_id, number_of_reactants,number_of_products,reactant_1,  reactant_2,  product_1,   product_2,  rate,         dG,          dG_barrier,  is_redox,
                            template,mapped_rxn,rxn)
                             values 
                            ({tmp_row[0]},{tmp_row[1]},       {tmp_row[2]},      {tmp_row[3]},{tmp_row[4]},{tmp_row[5]},{tmp_row[6]},{tmp_row[7]},{tmp_row[8]},{tmp_row[9]},{tmp_row[10]},
                            {template},{mapped_rxn},{tmp_rxn})
                            """
                filtered_rxn_cur.execute(filtered_sql_str)

    filtered_rxn_con.commit()

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
        well_define = int(row["well_define"])
        mol_charge = int(row["mol_charge"])
        star_atom_num=int(row["star_atom_num"])
        if well_define==1 and mol_charge==0 and star_atom_num==0:
            id_smiles_dict[mol_id] = smiles
            smiles_id_dict[smiles] = mol_id


    rn_con = sqlite3.connect(rn_db_path)
    rn_cur = rn_con.cursor()

    sql_limit = 100000
    rn_cur.execute(
        f"select  reaction_id, number_of_reactants, number_of_products, reactant_1, reactant_2, product_1, product_2 from "
        f"reactions where reaction_id%10000=0 limit {sql_limit}")

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
    """
    https://github.com/neo-chem-synth-wave/atom-to-atom-mapping
    """
    #add columns mapped_reaction_smiles,mapped_reaction_template_smarts,is_confident
    os.system(f"""python atom-to-atom-mapping/scripts/map_reaction_smiles_using_local_mapper.py \
     --batch_size 32 --input_csv_file_path {rxn_smarts_file}  \
     --reaction_smiles_column_name  rxn_smarts \
     --output_csv_file_path  {mapped_rxn_smarts_output_file}""")

def eda_mapped_rxn_smarts(mapped_rxn_smarts_file: str):
    df=pd.read_csv(mapped_rxn_smarts_file)

    output_dir=f"{data_dir}tmp_rxn"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    height_mol = 300
    width_mol = 300
    template_dict={}
    for i,( _,row) in enumerate(tqdm(df.iterrows(),total=df.shape[0])):
        rxn_smarts=row["rxn_smarts"]
        mapped_reaction_smiles=row["mapped_reaction_smiles"]
        mapped_reaction_template_smarts=row["mapped_reaction_template_smarts"]
        is_confident=row["is_confident"]
        rxn_id = row["reaction_id"]
        # draw_reaction(rxn_smarts,f"{output_dir}/{i}_{is_confident}.png")
        # draw_reaction(mapped_reaction_smiles,f"{output_dir}/{i}_{is_confident}_mapped.png")
        # draw_reaction(mapped_reaction_template_smarts,f"{output_dir}/{i}_{is_confident}_template.png")

        if is_confident==False:
            k="false"
        else:
            k=mapped_reaction_template_smarts
        if k not in template_dict:
            template_dict[k]=0
        template_dict[k]+=1
        if template_dict[k] < 4 :
            if_draw=True
        else:
            if_draw = False


        if if_draw:
            pil_img_list=[]
            # rxn = AllChem.ReactionFromSmarts(rxn_smarts, useSmiles=True)
            # img = ReactionToImage(rxn)
            img=draw_reaction(rxn_smarts,save=False)
            img.resize((width_mol*5, height_mol))
            pil_img_list.append(img)
            rxn = AllChem.ReactionFromSmarts(mapped_reaction_smiles, useSmiles=True)
            img = ReactionToImage(rxn)
            img.resize((width_mol*5, height_mol))
            pil_img_list.append(img)
            try:
                rxn = AllChem.ReactionFromSmarts(mapped_reaction_template_smarts, useSmiles=False)
                img = ReactionToImage(rxn)
                img.resize((width_mol*5, height_mol))
                pil_img_list.append(img)
            except:
                pass
            mode = pil_img_list[0].mode
            # 创建一个空白画布，用于拼接图片
            result = Image.new(mode, (width_mol * 5, height_mol*len(pil_img_list)), color=(255, 255, 255))  #

            # 在画布上拼接图片
            for img_idx, img in enumerate(pil_img_list):
                if img is None:
                    continue
                result.paste(img, (0, height_mol * img_idx, ))

            result.save(f"{output_dir}/{rxn_id}_{is_confident}.png")

def apply_MechFinder(mapped_rxn_smarts_file: str, mech_output_file: str):
    #from rxnmapper import RXNMapper  这句话会导致MechFinder无法打印报错
    from MechFinder import MechFinder
    finder = MechFinder(collection_dir='MechFinder/collections')
    df=pd.read_csv(mapped_rxn_smarts_file)

    output_dir=f"{data_dir}tmp_rxn"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    result_info_dict={}
    with open(mech_output_file, "w") as fp_out:
        print("reaction_id,rxn_str,updated_reaction,LRT,MT_class,electron_path",file=fp_out)
        # df = df[df['is_confident'] == "True"]
        for i,( _,row) in enumerate(tqdm(df.iterrows(),total=df.shape[0])):
            if i%1000==0:
                print(result_info_dict)
            if i>4:
                break
            rxn_id=row["reaction_id"]
            mapped_rxn = row["mapped_reaction_smiles"]
            unmapped_rxn = row["rxn_smarts"]
            #mapping_confidence=row["confidence"]

            try:
                updated_reaction, LRT, MT_class, electron_path = finder.get_electron_path(mapped_rxn)
            except Exception as e :
                if "except" not in result_info_dict:
                    result_info_dict["except"]=0
                result_info_dict["except"]+=1
                MT_class="except"
            else:
                if MT_class not in result_info_dict:
                    result_info_dict[MT_class] = 0
                result_info_dict[MT_class] += 1
                if not isinstance(finder.check_exception(MT_class), str):
                    if MT_class!="mechanism not in collection":
                        print(f"{rxn_id},{mapped_rxn},{updated_reaction},{LRT},{MT_class},{electron_path}",file=fp_out)
            #print(MT_class)
            draw_reaction(mapped_rxn, f"{output_dir}/{rxn_id}.png")
            print(f"{rxn_id},{mapped_rxn},{unmapped_rxn},{MT_class}")

def count_elements(mol):
    """
    统计分子中各类元素的数量。

    参数:
    mol (rdkit.Chem.Mol): RDKit 分子对象

    返回:
    dict: 元素及其对应的数量
    """
    element_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in element_counts:
            element_counts[symbol] += 1
        else:
            element_counts[symbol] = 1
    return element_counts

def draw_molecule(smiles, filename="molecule.png"):
    """
    从 SMILES 字符串绘制分子结构图并保存为图片文件。

    参数:
    smiles (str): 分子的 SMILES 表示
    filename (str): 保存图片的文件名，默认为 "molecule.png"
    """
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    img.save(filename)

def draw_reaction(rxn_smarts, filename="reaction.png",save=True):
    # rxn = AllChem.ReactionFromSmarts(rxn_smarts, useSmiles=True)
    # img = ReactionToImage(rxn)
    # img.save(filename)

    reactant_str, product_str = rxn_smarts.split(">>")

    reactant_list = reactant_str.split(".")
    reactant_1=reactant_list[0]
    if len(reactant_list)==2:
        reactant_2=reactant_list[1]
    else:
        reactant_2=""
    product_list = product_str.split(".")
    product_1=product_list[0]
    if len(product_list)==2:
        product_2=product_list[1]
    else:
        product_2=""
    smiles_list = [reactant_1, reactant_2,"", product_1, product_2]
    mode=None
    pil_img_list=[]
    for idx, smiles in enumerate(smiles_list):
        if smiles == "":
            img=None
        else:
            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol)
            mode=img.mode
        pil_img_list.append(img)

    height_mol = 300
    width_mol = 300
    # 创建一个空白画布，用于拼接图片
    result = Image.new(mode, (width_mol*5, height_mol), color=(255, 255, 255))  #

    # 在画布上拼接图片
    for img_idx, img in enumerate(pil_img_list):
        if img is None:
            continue
        result.paste(img, (width_mol*img_idx, 0))

    # 保存拼接后的图片
    if save:
        result.save(filename)
    else:
        return result
def find_reaction(smi_csv_path: str,rn_db_path: str):

    output_dir=f"{data_dir}tmp_rxn"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    id_smiles_dict = {}
    smiles_id_dict = {}
    smi_df = pd.read_csv(smi_csv_path)
    for i, row in tqdm(smi_df.iterrows(), total=smi_df.shape[0]):
        mol_id = int(row["idx"])
        smiles = row["smiles"]
        id_smiles_dict[mol_id] = smiles
        smiles_id_dict[smiles] = mol_id


    rn_con = sqlite3.connect(rn_db_path)
    rn_cur = rn_con.cursor()

    rn_cur.execute(
        f"select  reaction_id, number_of_reactants, number_of_products, reactant_1, reactant_2, product_1, product_2 from "
        f"reactions where  number_of_reactants=2 and number_of_products=1 and  product_1=13203 ") #and (reactant_1=13590 or reactant_2=13590)

    for row in tqdm(rn_cur):
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

        rxn_id_rxn_str=f"{reactant_1}.{reactant_2}>>{product_1}.{product_2}"
        print(f"{reaction_id},{rxn_smarts},{rxn_id_rxn_str}" )
        draw_reaction(rxn_smarts, f"{output_dir}/{reaction_id}_{reactant_1}.{reactant_2}>>{product_1}.{product_2}_{rxn_id_rxn_str}.png")

def find_mol(smiles_csv_file:str):
    from rdkit.Chem import Descriptors
    output_dir=f"{data_dir}tmp"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

    df=pd.read_csv(smiles_csv_file)
    for i,row in tqdm(df.iterrows(),total=df.shape[0]):
        smiles=row["smiles"]
        well_define=row["well_define"]
        mol_id=row["idx"]
        mol=Chem.MolFromSmiles(smiles)
        if mol:
            count_elements_dict=count_elements(mol)
            ring_count = Descriptors.RingCount(mol)
            # 获取分子中的所有环
            sssr = mol.GetRingInfo().AtomRings()
            # 检查是否有五元环
            has_five_membered_ring = any(len(ring) == 5 for ring in sssr)
            # 检查分子中是否有自由基
            has_radical = any(atom.GetNumRadicalElectrons() > 0 for atom in mol.GetAtoms())
            # 确保分子的电荷被正确计算
            AllChem.ComputeGasteigerCharges(mol)

            # 计算分子的净电荷
            net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

            try:
                #if count_elements_dict["C"]==14 and count_elements_dict["O"]==2 and count_elements_dict["Li"]==2 and count_elements_dict["F"]==2:
                if net_charge==0 and well_define==1 and not has_radical and has_five_membered_ring and ring_count==1 and \
                        count_elements_dict["C"]==3 and count_elements_dict["O"]==3 and count_elements_dict["F"]==1  \
                        and "Li" not in count_elements_dict and "P" not in count_elements_dict:
                    print(f"{mol_id},{smiles},{well_define}")
                    draw_molecule(smiles,f"{output_dir}/{mol_id}.png")
            except:
                continue


if __name__ == "__main__":

    data_dir="/personal/Bohrium_task_hiprgen_rn/hiprgen_json2rn_output/libe_and_fmol_0911_all/"
    mol_entries_file=f"{data_dir}mol_entries.pickle"

    #find_mol(f"{data_dir}smiles.csv")
    #find_reaction(f"{data_dir}smiles.csv",f"/root/HiPRGen/data/libe_and_fmol_0911_all/rn.sqlite")

    #filter_mol_entries(mol_entries_file, f"{data_dir}smiles.csv")
    # trans_rxn_db2smarts(f"{data_dir}smiles.csv",
    #                     rn_db_path="/root/HiPRGen/data/libe_and_fmol_0911_all/rn.sqlite",
    #                     rxn_smarts_output_file=f"{data_dir}rxn_smarts.csv")
    #mapping_rxn(f"{data_dir}rxn_smarts.csv",f"{data_dir}rxn_smarts_mapped.csv")
    #apply_MechFinder(f"{data_dir}rxn_smarts_mapped.csv",f"{data_dir}rxn_smarts_mapped_mech.csv")
    #eda_mapped_rxn_smarts(f"{data_dir}rxn_smarts_mapped.csv")
    filter_rxn(f"{data_dir}smiles.csv",f"/root/HiPRGen/data/libe_and_fmol_0911_all/rn.sqlite",f"{data_dir}/rn_filtered.sqlite")