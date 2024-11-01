


#  放置到路径 /personal/Bohrium_task_hiprgen_rn
import json
import os
import shutil
import sqlite3
import pickle
from tqdm import tqdm



def split_rxn_db(original_rxn_db_path,split_num,output_dir):

    # 连接到 SQLite 数据库
    ori_conn = sqlite3.connect(original_rxn_db_path)
    ori_cursor = ori_conn.cursor()
    # 计算原表的行数
    ori_cursor.execute(f"SELECT COUNT(*) FROM reactions")
    total_rows = ori_cursor.fetchone()[0]
    print(f"Total number of rows in original table: {total_rows}")
    # 计算每个新表应该有多少行
    split_table_rows = total_rows // split_num
    print(f"Each split table will have {split_table_rows} rows")
    ori_offset = 0
    for split_id in range(split_num):
        if os.path.exists(f"{output_dir}/{split_id}_rn.sqlite"):
            print(f"{output_dir}/{split_id}_rn.sqlite exists")
            continue
        split_conn=sqlite3.connect(f"{output_dir}/{split_id}_rn.sqlite")
        split_cursor = split_conn.cursor()
        print(f"Creating split table {split_id}")
        split_cursor.execute("""
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
        split_conn.commit()

        if split_id==split_num-1:
            split_table_rows=total_rows-ori_offset
        orl_sql_str = f"""select * from reactions limit {split_table_rows} offset {ori_offset} ;"""
        print(orl_sql_str)
        ori_cursor.execute(orl_sql_str)
        for row in tqdm(ori_cursor):
            split_cursor.execute("INSERT INTO reactions VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", row)

        split_conn.close()
        ori_offset += split_table_rows

    ori_conn.close()
    print("done")




def submit_task(machine_num, input_dir, output_root_dir, image_name):




    for machine_id in range(machine_num):
        if machine_id in []:
            continue

        split_task_output_dir = f"{output_root_dir}/{machine_num}_{machine_id}"
        if os.path.exists(split_task_output_dir):
            print(f"split_task_output_dir exists")
            continue

        split_rxn_db_input_filename = f"rn_{machine_id}.sqlite"
        split_filtered_rxn_db_output_filename= f"rn_filtered_{machine_id}.sqlite"
        python_str=f" python /root/HiPRGen/MechFinder_rxn_filter_study/filter_rxn_mpi_start.py  {split_rxn_db_input_filename} {split_filtered_rxn_db_output_filename} "
        lbg_task_json_dict={
                "job_name": f"hiprgen_rn_filter_libe_fmol_{machine_num}_{machine_id}",
                "command":python_str,
                "platform": "ali",
                "disk_size": 200,
                "machine_type": "c128_m512_cpu",
                "image_name": image_name,
                "program_id": 14480,
                "input":input_dir,
                "result":split_task_output_dir,
        }
        json_params_file_path=f"{input_dir}/lbg_task_{machine_num}_{machine_id}.json"
        json.dump(lbg_task_json_dict, open(json_params_file_path, "w"), indent=2)
        print(json.dumps(lbg_task_json_dict, indent=2))
        shell_str=f"lbg job submit -i {json_params_file_path}"
        os.system(shell_str)
        print(shell_str)

def post_process(machine_num, output_root_dir, ):

    all_rn_conn = sqlite3.connect( f"{output_root_dir}/rn.sqlite")
    all_rn_cur = all_rn_conn.cursor()
    all_rn_cur.execute("""
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
        """ )
    all_rn_conn.commit()


    for machine_id in range(machine_num):
        split_task_output_dir = f"{output_root_dir}/{machine_num}_{machine_id}"
        split_filtered_rxn_db_output_filename = f"rn_filtered_{machine_id}.sqlite"
        split_conn = sqlite3.connect( f"{split_task_output_dir}/{split_filtered_rxn_db_output_filename}")
        split_cur = split_conn.cursor()

        split_sql_str = f"""select * from reactions;"""
        print(split_sql_str)
        split_cur.execute(split_sql_str)

        for row in tqdm(split_cur):
            all_rn_cur.execute("INSERT INTO reactions VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", row)
            all_rn_conn.commit()

        split_conn.close()
        all_rn_conn.close()


def main():
    machine_num = 30
    input_dir = "/personal/Bohrium_task_hiprgen_rn/hiprgen_json2rn_input/libe_and_fmol_0911_rn_filter/"
    if not os.path.exists(input_dir):
        os.makedirs(input_dir)
    output_root_dir = "/personal/Bohrium_task_hiprgen_rn/hiprgen_json2rn_output/libe_and_fmol_0911_rn_filter/"
    if not os.path.exists(output_root_dir):
        os.makedirs(output_root_dir)

    image_name = "registry.dp.tech/dptech/prod-17396/hiprgen:20241101"
    original_rxn_db_path = f"/root/HiPRGen/data/libe_and_fmol_0911_all/rn.sqlite"

    split_rxn_db(original_rxn_db_path, machine_num, input_dir)

    submit_task(machine_num, input_dir, output_root_dir, image_name)
    #post_process(machine_num, output_root_dir)


if __name__ == '__main__':
    main()