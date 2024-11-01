import subprocess
import os
import logging
import sys

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s %(module)s: %(lineno)d %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)



number_of_threads = str(int(int(os.popen("nproc").read().strip())*1.5))
logger.info(f"number_of_threads:{number_of_threads}")

rxn_db_path=sys.argv[1]#f"/root/HiPRGen/data/libe_and_fmol_0911_all/rn.sqlite"
filtered_rxn_db_path_path=sys.argv[2]#f"rn_filtered.sqlite"

subprocess.run(
    [
        'mpirun',
        #'--use-hwthread-cpus',
        "--oversubscribe",
        '-n',
        number_of_threads,
        'python',
        '/root/HiPRGen/MechFinder_rxn_filter_study/filter_rxn_mpi.py',
        f"/root/HiPRGen/data/libe_and_fmol_0911_all/mol_entries.pickle",
        rxn_db_path,
        filtered_rxn_db_path_path
    ]
)