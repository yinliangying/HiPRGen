# using dpdispatcher to perform multiple xtb calculations in a single node
import math
import os
import shutil
import random

import numpy as np
from ase.db.core import connect
from ase.io.gaussian import read_gaussian_out, write_gaussian_in
from dpdispatcher import Task, Submission, Machine, Resources


def local_gaussian(n_parallel_jobs, n_cpu_per_job, db_path, gaussian_key_line, cmd_line, **kwargs):
    cwd_ = os.getcwd()
    db_path = os.path.abspath(db_path)
    # prepare
    cooking_path = os.path.abspath('cooking')
    os.makedirs(cooking_path)
    mach_para = {
        'batch_type': "Shell",
        'context_type': "LazyLocalContext",
        'remote_root': '/root/test_dpdispatcher',
        'remote_profile': {},
        'local_root': cooking_path
    }

    resrc_para = {
        'number_node': 1,
        'cpu_per_node': n_cpu_per_job,
        'gpu_per_node': 0,
        'group_size': 1,
        'queue_name': "slurm",
        'envs': {
            'dummy': 'dummy/path',
        },
    }

    task_list = []
    id_name_list = []
    with connect(db_path) as db:
        for a_row in db.select():
            os.chdir(cooking_path)
            id_name = f'id_{a_row.real_id}_spin_{a_row.real_spin}_charge_{a_row.real_charge}'
            os.makedirs(id_name)
            os.chdir(id_name)
            id_name_list.append(id_name)
            an_atoms = a_row.toatoms()
            with open(file='gau.gjf', mode='w') as f:
                write_gaussian_in(fd=f,
                                  atoms=an_atoms,
                                  properties=[' '],
                                  method='',
                                  basis=gaussian_key_line,
                                  nprocshared=str(n_cpu_per_job),
                                  mem='10GB',  # 10GB by default
                                  mult=a_row.real_spin,
                                  charge=a_row.real_charge,
                                  chk='gau.chk',
                                  )
            with open('gau.gjf', "r+") as file:  # Open in read-write mode
                lines = file.readlines()  # Read all lines into a list
                del lines[-1]
                file.seek(0)  # Move the file pointer to the beginning
                file.writelines(lines)
                file.writelines(['Eps=18.5\n',
                                 'EpsInf=1.415\n',
                                 'HbondAcidity=0\n',
                                 'HbondBasicity=0.735\n',
                                 'SurfaceTensionAtInterface=20.2\n',
                                 'CarbonAromaticity=0\n',
                                 'ElectronegativeHalogenicity=0\n\n\n',
                                 ])
            a_task = Task(command=cmd_line,
                          task_work_path=f"{cooking_path}/{id_name}/",
                          forward_files=[f'{cooking_path}/{id_name}/*'],
                          backward_files=[]
                          )
            task_list.append(a_task)

    n_batches = math.ceil(len(task_list) / n_parallel_jobs)
    for i in range(n_batches):
        cursor = i * n_parallel_jobs
        a_task_list = task_list[cursor:cursor + n_parallel_jobs]
        submission = Submission(
            work_base=cooking_path,
            machine=Machine.load_from_dict(machine_dict=mach_para),
            resources=Resources.load_from_dict(resources_dict=resrc_para),
            task_list=a_task_list,
        )
        try:
            submission.run_submission()
        except:
            continue
    os.chdir(cwd_)


def remote_gaussian(n_parallel_machines, main_db_path, resrc_info, machine_info, handler_file_path):
    cwd_ = os.getcwd()
    # prepare
    abs_handler_file_path = os.path.abspath(handler_file_path)
    handler_file_name = os.path.basename(handler_file_path)
    cooking_path = os.path.abspath('cooking')
    if os.path.exists(cooking_path):
        print('Found previous workbase. It will be cleared.')
        shutil.rmtree(cooking_path)
    os.makedirs(cooking_path)
    task_list = []

    with connect(main_db_path) as main_db:
        total_n_mols = main_db.count()
        random_idx_list = list(range(total_n_mols))
        random.shuffle(random_idx_list)
        sub_n_mols = math.ceil(total_n_mols / n_parallel_machines)
        actual_machine_used = total_n_mols // sub_n_mols + 1
        for i in range(actual_machine_used):
            os.chdir(cooking_path)
            os.makedirs(f'{str(i)}')
            os.chdir(f'{str(i)}')
            shutil.copy(src=abs_handler_file_path, dst=handler_file_name)
            if i < actual_machine_used - 1:
                with connect('raw.db') as dump_db:
                    for row_idx, a_row in enumerate(main_db.select()):
                        if row_idx in random_idx_list[i * sub_n_mols: (i + 1) * sub_n_mols]:
                            dump_db.write(a_row)
            else:
                with connect('raw.db') as dump_db:
                    for row_idx, a_row in enumerate(main_db.select()):
                        if row_idx in random_idx_list[i * sub_n_mols: ]:
                            dump_db.write(a_row)
            # task
            a_task = Task(
                command=fr'unset SLURM_NTASKS && unset SLURM_JOB_NAME && python {handler_file_name} 2>&1 ',
                task_work_path=f'{str(i)}/',
                forward_files=[f'{cooking_path}/{str(i)}/*'],
                backward_files=[f'cooking']
            )
            task_list.append(a_task)
    os.chdir(cwd_)
    # submission
    submission = Submission(
        work_base=cooking_path,
        machine=Machine.load_from_dict(machine_dict=machine_info),
        resources=Resources.load_from_dict(resources_dict=resrc_info),
        task_list=task_list,
    )
    submission.run_submission()
    # # post-process
    # with connect(main_db_path) as old_db, connect('cooked.db') as new_db:
    #     print(f'The original db has {old_db.count()} molecules')
    #     for a_row in old_db.select():
    #         an_id = a_row.real_id
    #         for i in range(n_parallel_machines):
    #             opted_file = os.path.join(cooking_path, f'{str(i)}', 'cooking', f'id_{str(an_id)}', 'xtbopt.xyz')
    #             if os.path.exists(opted_file):
    #                 a_new_atoms = read(opted_file)
    #                 new_db.write(a_new_atoms, real_id=an_id, real_charge=a_row.real_charge, real_spin=a_row.real_spin)
    #     print(f'The optimized db has {new_db.count()} molecules')


