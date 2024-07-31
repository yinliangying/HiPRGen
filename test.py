import json
import os
machine_num=2

for machine_id in range(machine_num):
    python_str=f" python /root/HiPRGen/json_2_rn_pipeline.py  -m ab_initio -j new_libe/libe_and_fmol.json  -o  new_libe --machine_num {machine_num} --machine_id {machine_id}"
    json_dict={
            "job_name": f"hiprgen_json2rn_{machine_num}_{machine_id}",
            "command":python_str,
            "platform": "ali",
            "disk_size": 200,
            "machine_type": "c64_m512_cpu",
            "image_name": "registry.dp.tech/dptech/prod-17396/hiprgen:20240728",
            "program_id": 14480
    }
    json.dump(json_dict,open("hiprgen_json2rn_input/bohrium_task_json2rn.json","w"),indent=2)
    print(json_dict)
    shell_str=f"lbg job submit -i hiprgen_json2rn_input/bohrium_task_json2rn.json  -p hiprgen_json2rn_input  -r hiprgen_json2rn_output/{machine_num}_{machine_id}"
    os.system(shell_str)
    print(shell_str)