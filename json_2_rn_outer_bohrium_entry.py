#  放置到路径 /personal/Bohrium_task_hiprgen_rn
import json
import os
import shutil
import sqlite3
import pickle


machine_num=10
json_input_file_name= "libe_and_fmol_0911.json"
json_input_file_dir= "hiprgen_json2rn_input/libe_and_fmol_0911/"
output_dir_prefix="hiprgen_json2rn_output/libe_and_fmol_0911_"
image_name="registry.dp.tech/dptech/prod-17396/hiprgen:202409142"

for machine_id in range(machine_num):
    python_str=f" python /root/HiPRGen/json_2_rn_bohirum_entry.py  -m ab_initio -j {json_input_file_name}  -o  ./ --machine_num {machine_num} --machine_id {machine_id}"
    json_dict={
            "job_name": f"hiprgen_json2rn_libe_fmol_{machine_num}_{machine_id}",
            "command":python_str,
            "platform": "ali",
            "disk_size": 200,
            "machine_type": "c64_m512_cpu",
            "image_name": image_name,
            "program_id": 14480,
            "input":json_input_file_dir,
            "result":f"{output_dir_prefix}{machine_num}_{machine_id}",
    }
    json_params_file_path=f"{json_input_file_dir}/lbg_task_{machine_num}_{machine_id}.json"
    json.dump(json_dict,open(json_params_file_path,"w"),indent=2)
    print(json.dumps(json_dict,indent=2))
    shell_str=f"lbg job submit -i {json_params_file_path}"
    os.system(shell_str)
    print(shell_str)

exit()

#全部运行完后 把各个子表合并到all中
output_dir_all=f"{output_dir_prefix}all"
if not os.path.exists(output_dir_all):
    os.mkdir(output_dir_all)

shutil.copy(f"{output_dir_prefix}{machine_num}_0/mol_entries.pickle",f"{output_dir_all}/mol_entries.pickle")
all_number_of_species=len(pickle.load(open(f"{output_dir_all}/mol_entries.pickle","rb")))

all_rn_con = sqlite3.connect( f"{output_dir_all}/rn.sqlite")
all_rn_cur = all_rn_con.cursor()
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
            is_redox            INTEGER NOT NULL
    );
""")
all_rn_con.commit()

limit=100000
all_reaction_id=0
for machine_id in range(machine_num):

    rn_con = sqlite3.connect( f"{output_dir_prefix}{machine_num}_{machine_id}/rn.sqlite")
    rn_cur = rn_con.cursor()

    machine_number_of_reactions = 0
    machine_sql_str=f"""select number_of_reactions  from  metadata;"""
    print(machine_sql_str)
    rn_cur.execute(machine_sql_str)
    for (tmp_number_of_reactions,) in rn_cur:
        machine_number_of_reactions=tmp_number_of_reactions
        break
    print(f" machine_id:{machine_id},machine_number_of_reactions:{machine_number_of_reactions}")


    for offset in range(0,machine_number_of_reactions,limit):
        machine_sql_str=f"""select 
         reaction_id ,number_of_reactants ,number_of_products , reactant_1  ,reactant_2 ,product_1 ,
         product_2  ,rate  , dG   ,dG_barrier  ,is_redox 
         from reactions limit {limit} offset {offset} ;
        """
        print(machine_sql_str)
        rn_cur.execute(machine_sql_str)

        for i,(_, number_of_reactants, number_of_products, reactant_1, reactant_2, product_1,
             product_2, rate, dG, dG_barrier, is_redox) in enumerate(rn_cur):


            all_sql_str=f"""
            insert into reactions 
            (reaction_id      ,  number_of_reactants ,  number_of_products ,  reactant_1  , reactant_2 , product_1 ,  product_2  , rate  , dG   ,dG_barrier  ,is_redox )
             values 
            ({all_reaction_id}, {number_of_reactants}, {number_of_products}, {reactant_1}, {reactant_2}, {product_1},{product_2}, {rate}, {dG}, {dG_barrier}, {is_redox})
            """
            if i==0:
                print(machine_id,all_sql_str)
            all_rn_cur.execute(all_sql_str)
            all_reaction_id+=1

        all_rn_con.commit()

all_sql_str="""
    CREATE TABLE metadata (
            number_of_species   INTEGER NOT NULL,
            number_of_reactions INTEGER NOT NULL
    );
"""
all_rn_cur.execute(all_sql_str)
all_sql_str=f"""
    insert into metadata (number_of_species,number_of_reactions) values ({all_number_of_species},{all_reaction_id})
"""
all_rn_cur.execute(all_sql_str)
all_rn_con.commit()