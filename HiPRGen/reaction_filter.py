from mpi4py import MPI
from itertools import permutations, product
from HiPRGen.report_generator import ReportGenerator
from time import time
from HiPRGen.logging import llog_message
import sqlite3
from enum import Enum
from math import floor

from HiPRGen.reaction_questions import (
    run_decision_tree
)

"""
Phases 3 & 4 run in paralell using MPI

Phase 3: reaction gen and filtering
input: a bucket labeled by atom count
output: a list of reactions from that bucket
description: Loop through all possible reactions in the bucket and apply the decision tree. This will run in parallel over each bucket.

Phase 4: collating and indexing
input: all the outputs of phase 3 as they are generated
output: reaction network database
description: the worker processes from phase 3 are sending their reactions to this phase and it is writing them to DB as it gets them. We can ensure that duplicates don't get generated in phase 3 which means we don't need extra index tables on the db.

the code in this file is designed to run on a compute cluster using MPI.
"""

def table_exists(connection, table_name):
    cursor = connection.cursor()
    query = """
    SELECT name FROM sqlite_master 
    WHERE type='table' AND name=?;
    """
    cursor.execute(query, (table_name,))
    result = cursor.fetchone()
    return result is not None

create_metadata_table = """
    CREATE TABLE metadata (
            number_of_species   INTEGER NOT NULL,
            number_of_reactions INTEGER NOT NULL
    );
"""

insert_metadata = """
    INSERT INTO metadata VALUES (?, ?)
"""

# it is important that reaction_id is the primary key
# otherwise the network loader will be extremely slow.
create_reactions_table = """
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
"""


insert_reaction = """
    INSERT INTO reactions VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
"""

get_complex_group_sql = """
    SELECT * FROM complexes WHERE composition_id=? AND group_id=?
"""


# TODO: structure these global variables better
DISPATCHER_RANK = 0

# message tags

# sent by workers to the dispatcher once they have finished initializing
# only sent once
INITIALIZATION_FINISHED = 0

# sent by workers to the dispatcher to request a new table
SEND_ME_A_WORK_BATCH = 1

# sent by dispatcher to workers when delivering a new table
HERE_IS_A_WORK_BATCH = 2

# sent by workers to the dispatcher when reaction passes db decision tree
NEW_REACTION_DB = 3

# sent by workers to the dispatcher when reaction passes logging decision tree
NEW_REACTION_LOGGING = 4

class WorkerState(Enum):
    INITIALIZING = 0
    RUNNING = 1
    FINISHED = 2

def dispatcher(
        mol_entries,
        dispatcher_payload
):

    comm = MPI.COMM_WORLD
    work_batch_list = []
    bucket_con = sqlite3.connect(dispatcher_payload.bucket_db_file)
    bucket_cur = bucket_con.cursor()
    size_cur = bucket_con.cursor()

    #整体bucket的结构是 bucket的id是composition_id,原理上桶内所有元素要两两组合，但考虑计算效率进行分组，count是桶内分组数
    res = bucket_cur.execute("SELECT * FROM group_counts")
    tmp_i=0


    for (composition_id, count) in res:
        for (i,j) in product(range(count), repeat=2): #桶内分组两两组合
            if (dispatcher_payload.machine_id is not None)   and (dispatcher_payload.machine_num is not  None):
                if tmp_i%dispatcher_payload.machine_num==dispatcher_payload.machine_id:
                    work_batch_list.append(
                        (composition_id, i, j))
            else:
                work_batch_list.append(
                    (composition_id, i, j))
            tmp_i+=1

    composition_names = {}
    res = bucket_cur.execute("SELECT * FROM compositions")
    for (composition_id, composition) in res:
        composition_names[composition_id] = composition

    log_message("creating reaction network db")
    rn_con = sqlite3.connect(dispatcher_payload.reaction_network_db_file)
    rn_cur = rn_con.cursor()
    rn_cur.execute(create_metadata_table)

    if not table_exists(rn_con, "reactions"):
        #从头建立数据模式
        log_message("creating reaction network db")
        rn_cur.execute(create_reactions_table)
    else:#追加数据模式
        pass
    rn_con.commit()

    log_message("initializing report generator")

    # since MPI processes spin lock, we don't want to have the dispathcer
    # spend a bunch of time generating molecule pictures
    report_generator = ReportGenerator(
        mol_entries,
        dispatcher_payload.report_file,
        rebuild_mol_pictures=False
    )

    worker_states = {}

    worker_ranks = [i for i in range(comm.Get_size()) if i != DISPATCHER_RANK]

    for i in worker_ranks:
        worker_states[i] = WorkerState.INITIALIZING

    for i in worker_states:
        # block, waiting for workers to initialize
        comm.recv(source=i, tag=INITIALIZATION_FINISHED)
        worker_states[i] = WorkerState.RUNNING

    log_message("all workers running")

    res = rn_cur.execute("select count(*) from reactions")
    for row in res:
        reaction_index = row[0]
    log_message("old reaction index:", reaction_index)

    log_message("handling requests")

    batches_left_at_last_checkpoint = len(work_batch_list)
    last_checkpoint_time = floor(time())
    while True:
        if WorkerState.RUNNING not in worker_states.values():
            break

        current_time = floor(time())
        time_diff = current_time - last_checkpoint_time
        if ( current_time % dispatcher_payload.checkpoint_interval == 0 and
             time_diff > 0):
            batches_left_at_current_checkpoint = len(work_batch_list)
            batch_count_diff = (
                batches_left_at_last_checkpoint -
                batches_left_at_current_checkpoint)

            batch_consumption_rate = batch_count_diff / time_diff * 60

            log_message("batches remaining:", batches_left_at_current_checkpoint)
            log_message("batch consumption rate:",
                        batch_consumption_rate,
                        "batches per minute")


            batches_left_at_last_checkpoint = batches_left_at_current_checkpoint
            last_checkpoint_time = current_time


        status = MPI.Status()
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()
        rank = status.Get_source()

        if tag == SEND_ME_A_WORK_BATCH:
            if len(work_batch_list) == 0:
                comm.send(None, dest=rank, tag=HERE_IS_A_WORK_BATCH)
                worker_states[rank] = WorkerState.FINISHED
            else:
                # pop removes and returns the last item in the list
                work_batch = work_batch_list.pop()
                comm.send(work_batch, dest=rank, tag=HERE_IS_A_WORK_BATCH)
                composition_id, group_id_0, group_id_1 = work_batch
                # log_message(
                #     "dispatched",
                #     composition_names[composition_id],
                #     ": group ids:",
                #     group_id_0, group_id_1
                # )


        elif tag == NEW_REACTION_DB:
            reaction = data
            rn_cur.execute(
                insert_reaction,
                (reaction_index,
                 reaction['number_of_reactants'],
                 reaction['number_of_products'],
                 reaction['reactants'][0],
                 reaction['reactants'][1],
                 reaction['products'][0],
                 reaction['products'][1],
                 reaction['rate'],
                 reaction['dG'],
                 reaction['dG_barrier'],
                 reaction['is_redox']
                 ))

            reaction_index += 1
            if reaction_index % dispatcher_payload.commit_frequency == 0:
                rn_con.commit()


        elif tag == NEW_REACTION_LOGGING:

            reaction = data[0]
            decision_path = data[1]

            report_generator.emit_verbatim(decision_path)
            report_generator.emit_reaction(reaction)
            report_generator.emit_bond_breakage(reaction)
            report_generator.emit_newline()



    log_message("finalzing database and generation report")
    rn_cur.execute(
        insert_metadata,
        (len(mol_entries),
         reaction_index)
    )


    report_generator.finished()
    rn_con.commit()
    bucket_con.close()
    rn_con.close()


def worker(
        mol_entries,
        worker_payload
):
    print("预加载cpp交互数据")
    from HiPRGen.fragment_matching_found_cpp import create_molecule_entry
    for i in range(len(mol_entries)):
        mol_entry_ctype = create_molecule_entry(mol_entries, i)
        mol_entries[i].mol_entry_ctype = mol_entry_ctype


    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    con = sqlite3.connect(worker_payload.bucket_db_file)
    cur = con.cursor()

    comm.send(None, dest=DISPATCHER_RANK, tag=INITIALIZATION_FINISHED)

    total_time = 0
    comm_time = 0
    sql_time = 0
    filter_time = 0
    batch_times = 0
    while True:
        if rank == 1:
            batch_times += 1
            # print(
            #     f"total_time:{total_time} comm_time:{comm_time}  sql_time:{sql_time}  filter_time:{filter_time}   batch_times:{batch_times}   ")

        s_time=time()

        t_time = time()
        comm.send(None, dest=DISPATCHER_RANK, tag=SEND_ME_A_WORK_BATCH)
        work_batch = comm.recv(source=DISPATCHER_RANK, tag=HERE_IS_A_WORK_BATCH)
        comm_time += time() - t_time

        if work_batch is None:
            break

        composition_id, group_id_0, group_id_1 = work_batch

        if group_id_0 == group_id_1:

            t_time = time()
            res = cur.execute(
                get_complex_group_sql,
                (composition_id, group_id_0))
            sql_time += time() - t_time

            bucket = []
            for row in res:
                bucket.append((row[0], row[1]))

            iterator = permutations(bucket, r=2)

        else:
            t_time = time()
            res_0 = cur.execute(
                get_complex_group_sql,
                (composition_id, group_id_0))
            sql_time += time() - t_time

            bucket_0 = []
            for row in res_0:
                bucket_0.append((row[0], row[1]))

            t_time = time()
            res_1 = cur.execute(
                get_complex_group_sql,
                (composition_id, group_id_1))
            sql_time += time() - t_time

            bucket_1 = []
            for row in res_1:
                bucket_1.append((row[0], row[1]))

            iterator = product(bucket_0, bucket_1)

        for (reactants, products) in iterator:
            reaction = {
                'reactants': reactants,
                'products': products,
                'number_of_reactants': len([i for i in reactants if i != -1]),
                'number_of_products': len([i for i in products if i != -1])}

            t1_time = time()

            decision_pathway = []
            if run_decision_tree(reaction,
                                 mol_entries,
                                 worker_payload.params,
                                 worker_payload.reaction_decision_tree,
                                 decision_pathway
                                 ):

                t_time = time()
                comm.send(
                    reaction,
                    dest=DISPATCHER_RANK,
                    tag=NEW_REACTION_DB)
                comm_time += time() - t_time

            # if "dG" in reaction:
            #     if "fragment_matching_found" in reaction:
            #         print(f"XXX\t{reaction['dG']}\t{reaction['fragment_matching_found']}")
            if run_decision_tree(reaction,
                                 mol_entries,
                                 worker_payload.params,
                                 worker_payload.logging_decision_tree):

                t_time = time()
                comm.send(
                    (reaction,
                     '\n'.join([str(f) for f in decision_pathway])
                     ),

                    dest=DISPATCHER_RANK,
                    tag=NEW_REACTION_LOGGING)
                comm_time += time() - t_time

            filter_time += time() - t1_time

        total_time += time() - s_time
