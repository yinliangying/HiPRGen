import time

import pandas as pd
import sqlite3
import os
from tqdm import tqdm
from localmapper import localmapper
import logging
import traceback
import sys
import pickle
from mpi4py import MPI
from enum import Enum
from research import filter_mol_entries

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s: %(lineno)d %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

class DispatcherWorkerProcess():
    def __init__(self,mol_entries_path, rxn_db_path, filtered_rxn_db_path_path  ):
        self.cnt=0
        self.cmt_freq=100
        self.mapping_batch_size=1000
        self.filtered_rxn_db_path_path=filtered_rxn_db_path_path
        self.mol_entries_path=mol_entries_path
        self.rxn_db_path=rxn_db_path
    def init_dispatcher(self):
        # prepare filtered_rxn_db
        if os.path.exists(self.filtered_rxn_db_path_path):
            os.remove(self.filtered_rxn_db_path_path)
        self.filtered_rxn_con = sqlite3.connect(self.filtered_rxn_db_path_path)
        self.filtered_rxn_cur = self.filtered_rxn_con.cursor()
        self.filtered_rxn_cur.execute("""
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
        self.filtered_rxn_con.commit()


    def data_gen(self):
        smi_csv_path = "smiles.csv"
        filter_mol_entries(self.mol_entries_path, smi_csv_path)

        # Load smiles
        id_smiles_dict = {}
        smiles_id_dict = {}
        smi_df = pd.read_csv(smi_csv_path)
        for i, row in tqdm(smi_df.iterrows(), total=smi_df.shape[0]):
            mol_id = int(row["idx"])
            smiles = row["smiles"]
            well_define = int(row["well_define"])
            mol_charge = int(row["mol_charge"])
            star_atom_num = int(row["star_atom_num"])
            if True:  # if well_define==1 and mol_charge==0 and star_atom_num==0:
                id_smiles_dict[mol_id] = smiles
                smiles_id_dict[smiles] = mol_id
        logger.info(f"smiles_id_dict length:{len(smiles_id_dict)}")

        rn_con = sqlite3.connect(self.rxn_db_path)
        rn_cur = rn_con.cursor()

        rn_cur.execute(
            f"select  reaction_id, number_of_reactants, number_of_products, reactant_1, reactant_2, product_1, product_2, "
            f"rate,dG,dG_barrier,is_redox from reactions")


        tmp_mapping_row_list = []
        tmp_mapping_rxn_list = []

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
            rxn_smiles = reactant_1_smiles
            if number_of_reactants == 2:
                rxn_smiles += "." + reactant_2_smiles
            rxn_smiles += f">>{product_1_smiles}"
            if number_of_products == 2:
                rxn_smiles += "." + product_2_smiles

            tmp_mapping_rxn_list.append(rxn_smiles)
            tmp_mapping_row_list.append(row)

            if len(tmp_mapping_row_list) == self.mapping_batch_size:
                input_batch = (tmp_mapping_row_list, tmp_mapping_rxn_list)
                yield input_batch
                tmp_mapping_row_list = []
                tmp_mapping_rxn_list=[]


        input_batch = (tmp_mapping_row_list, tmp_mapping_rxn_list)
        yield input_batch
    def init_worker(self):
        self.mapper = localmapper("cpu")

    def worker_run(self,input_batch):
        tmp_mapping_row_list, tmp_mapping_rxn_list = input_batch
        try:
            result=self.mapper.get_atom_map(tmp_mapping_rxn_list, return_dict=True)
        except:
            logger.error(f"worker_run error ")
            logger.error(str(input_batch))
            logger.error(traceback.format_exc())
            result=None
        data=(tmp_mapping_row_list, tmp_mapping_rxn_list,result)
        return data

    def post_process(self,data):
        tmp_mapping_row_list, tmp_mapping_rxn_list, tmp_result_list=data
        if tmp_result_list==None:
            return
        for tmp_row,tmp_rxn, tmp_result in zip(tmp_mapping_row_list,tmp_mapping_rxn_list, tmp_result_list):
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
                    sql_data = list(tmp_row) + [template, mapped_rxn, tmp_rxn]
                    try:
                        self.filtered_rxn_cur.execute(filtered_sql_str, sql_data)
                    except:

                        logger.error(traceback.format_exc())
                        logger.error(str(sql_data))
                        continue
        self.cnt+=1
        if self.cnt%self.cmt_freq==0:
            self.filtered_rxn_con.commit()
            logger.info(f"commit {self.cnt}")


class MPItask():
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
    def __init__(self,dispatcher_worker_process):
        # TODO: structure these global variables better
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.dispatcher_worker_process=dispatcher_worker_process
        if self.rank == self.DISPATCHER_RANK:
            self.dispatcher()
        else:
            self.worker()
    def dispatcher(self,):
        gen=self.dispatcher_worker_process.data_gen()
        self.dispatcher_worker_process.init_dispatcher()

        worker_states = {}
        worker_ranks = [i for i in range(self.comm.Get_size()) if i != self.DISPATCHER_RANK]
        for i in worker_ranks:
            worker_states[i] = self.WorkerState.INITIALIZING
        for i in worker_states:
            # block, waiting for workers to initialize
            self.comm.recv(source=i, tag=self.INITIALIZATION_FINISHED)
            worker_states[i] = self.WorkerState.RUNNING


        while True:
            if self.WorkerState.RUNNING not in worker_states.values():
                time.sleep(1)
                break

            status = MPI.Status()
            data = self.comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            tag = status.Get_tag()
            worker_rank = status.Get_source()

            if tag == self.SEND_ME_A_WORK_BATCH:
                try:
                    input_batch = next(gen)
                except StopIteration:
                    # 生成器结束，处理结束逻辑
                    logger.info("生成器已结束，没有更多数据可生成")
                    self.comm.send(None, dest=worker_rank, tag=self.HERE_IS_A_WORK_BATCH)
                    worker_states[worker_rank] = self.WorkerState.FINISHED

                else:
                    self.comm.send(input_batch, dest=worker_rank, tag=self.HERE_IS_A_WORK_BATCH)
            elif tag == self.NEW_REACTION_DB:
                self.dispatcher_worker_process.post_process(data)


    def worker(self):

        self.comm.send(None, dest=self.DISPATCHER_RANK, tag=self.INITIALIZATION_FINISHED)
        self.dispatcher_worker_process.init_worker()
        while True:
            self.comm.send(None, dest=self.DISPATCHER_RANK, tag=self.SEND_ME_A_WORK_BATCH)
            input_batch = self.comm.recv(source=self.DISPATCHER_RANK, tag=self.HERE_IS_A_WORK_BATCH)
            if input_batch is None:
                break
            self.comm.send(
                self.dispatcher_worker_process.worker_run(input_batch),
                dest=self.DISPATCHER_RANK,
                tag=self.NEW_REACTION_DB)

#
# def test_yield():
#     print("init")
#     for i in range(1000000):
#         yield i
#     yield "end"


def main():
    mol_entries_path = sys.argv[1]
    rxn_db_path = sys.argv[2]
    filtered_rxn_db_path_path = sys.argv[3]
    dispatcher_worker_process=DispatcherWorkerProcess(mol_entries_path,rxn_db_path,filtered_rxn_db_path_path)
    MPItask( dispatcher_worker_process)
    # gen=dispatcher_worker_process.data_gen()
    # while True:
    #     #input_batch=test_yield()
    #     input_batch = next(gen)
    #     print(input_batch)



if __name__ == "__main__":

    main()


