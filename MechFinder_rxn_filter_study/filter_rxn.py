
import pandas as pd
import sqlite3
import os
from tqdm import tqdm
from localmapper import localmapper
import logging
import traceback

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s: %(lineno)d %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)




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
        if True:#if well_define==1 and mol_charge==0 and star_atom_num==0:
            id_smiles_dict[mol_id] = smiles
            smiles_id_dict[smiles] = mol_id
    logger.info(f"smiles_id_dict length:{len(smiles_id_dict)}")
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
    mapping_batch_size=1000
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
            logger.info(f"mapping times:{mapping_times}")
            tmp_result_list = mapper.get_atom_map(tmp_mapping_rxn_list, return_dict=True)
            logger.info("mapping finished")
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
                        (?,?,?,?,?,?,?,?,?,?,?,
                        ?,?,?)
                        """
                        try:
                            sql_data = list(tmp_row) + [template, mapped_rxn, tmp_rxn]
                            filtered_rxn_cur.execute(filtered_sql_str, sql_data)
                        except:

                            logger.error(traceback.format_exc())
                            logger.error(str(sql_data))
                            continue
            tmp_mapping_rxn_list=[]
            tmp_mapping_row_list=[]

            if mapping_times % commit_freq == 0:
                logger.info(f"commit at {mapping_times}")
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
                (?,?,?,?,?,?,?,?,?,?,?,
                ?,?,?)
                """
                try:
                    sql_data=list(tmp_row)+[template,mapped_rxn,tmp_rxn]
                    filtered_rxn_cur.execute(filtered_sql_str, sql_data)
                except:

                    logger.error(traceback.format_exc())
                    logger.error(str(sql_data))
                    continue

    filtered_rxn_con.commit()




if __name__ == "__main__":

    data_dir="/personal/Bohrium_task_hiprgen_rn/hiprgen_json2rn_output/libe_and_fmol_0911_all/"
    mol_entries_file=f"{data_dir}mol_entries.pickle"

    filter_rxn(f"{data_dir}smiles.csv",f"/root/HiPRGen/data/libe_and_fmol_0911_all/rn.sqlite",f"{data_dir}/rn_filtered.sqlite")
