import os
import pickle
import sqlite3
import subprocess
import sys

import pandas
from HiPRGen.bucketing import bucket
from HiPRGen.constants import ROOM_TEMP, Terminal
from HiPRGen.initial_state import find_mol_entry_from_xyz_and_charge, find_mol_entry_by_entry_id
from HiPRGen.initial_state import insert_initial_state
from HiPRGen.mc_analysis import (
    reaction_tally_report,
    species_report,
    Pathfinding,
    SimulationReplayer,
    generate_pathway_report,
    sink_report,
    consumption_report,
    redox_report,
    coordination_report,
    decoordination_report
)
from HiPRGen.network_loader import NetworkLoader
from HiPRGen.reaction_filter_payloads import (
    DispatcherPayload,
    WorkerPayload
)
from HiPRGen.reaction_questions import (
    default_reaction_decision_tree,

)
from HiPRGen.report_generator import ReportGenerator
from HiPRGen.species_filter import species_filter
from HiPRGen.species_questions import (
    mg_species_decision_tree,
    li_species_decision_tree,
    positive_penalty,
    species_default_true
)
from monty.serialization import loadfn, dumpfn
from rdkit import Chem
from rdkit.Chem import AllChem

import logging

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] [%(name)s]  %(funcName)s  %(lineno)d- %(message)s')

logger = logging.getLogger(__name__)
def smile_to_inchi(smile):
    mol = Chem.MolFromSmiles(smile)
    inchi = Chem.MolToInchi(mol)
    return inchi


def smiles_to_xyz(smiles: str, output_file_path: str):
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)

    # Add hydrogen atoms
    mol = Chem.AddHs(mol)

    # Generate a 3D conformation
    AllChem.EmbedMolecule(mol)

    # Optimize the geometry
    AllChem.MMFFOptimizeMolecule(mol)

    # Get the atomic positions
    conf = mol.GetConformer()
    atoms = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]

    # Write the XYZ file
    xyz_content = f'{mol.GetNumAtoms()}\n\n'
    for atom, pos in zip(mol.GetAtoms(), atoms):
        symbol = atom.GetSymbol()
        x, y, z = pos
        xyz_content += f'{symbol} {x:10.6f} {y:10.6f} {z:10.6f}\n'
    with open(output_file_path, 'w') as f:
        f.write(xyz_content)


# Since HiPRGen uses an end-to-end testing approach rather than testing
# each individual function, we have decided to use the tests as
# documentation, by explaining every single line through the first test.


# The first thing you need to consider when using HiPRGen is how many
# worker threads do you want to run. HiPRGen can be run with a single
# thread or thousands distrubuted across several nodes. For reaction
# networks with between ~5000 and ~10000 species, we have found that the
# optimal number of worker threads is between 1000 and 2000. If you try
# and use more than that, the worker threads are going to spend lots of
# time waiting for the dispatcher to get through all of the reactions it
# is being sent, which slows everything down. Fixing this would require
# a more complex distrubuted system, but it hasn't been an issue even
# for very large reaction networks.


class bcolors:
    PASS = '\u001b[32;1m'
    FAIL = '\u001b[31;1m'
    ENDC = '\u001b[0m'


# HiPRGen is organized as a pipeline, where all the relevent data is
# stored in a sqlite database between phases. For this reason, during
# a single run of the full pipeline, it makes sense to store all the
# relevent files in a single directory. We have two test sets, a lithium
# set and a magnesium set. Since the lithium test set is older, we shall
# document that instead of the mg test set.

# if os.path.isdir('./scratch'):
#    subprocess.run(['rm', '-r', './scratch'])

# subprocess.run(['mkdir', './scratch'])


def gen_network(network_json_path,network_folder):

    if not os.path.exists(network_folder):
        os.system("mkdir -p %s" % (network_folder))
    # The initial input to the pipeline is a list of LIBE or MADEIRA
    # dataset entries. We provide two examples in the data foloder.

    database_entries = loadfn(network_json_path)

    # The first step of the HiPRGen pipeline is passing the input molecules
    # through the species decision tree to discard molecules.
    species_decision_tree = li_species_decision_tree

    # There is one non-local part of species filtering: we consider two
    # molecules to be equivalent if they have the same total charge,
    # composition, and covalent bonds, even if they have different metal
    # coordination, and we choose one such molecule in each "coordimer"
    # class using the coodimer weight function.  Since most of our logging
    # later on is defined in terms of a fixed molecule set, logging for
    # the species filtering phase is messy, so ignore the species_report
    # argument for now. The second argument is where we store a pickle of
    # the filtered molecule entries for use in later phases.
    # if os.path.exists(network_folder + '/mol_entries.pickle'):
    #     os.system("rm %s" % network_folder + '/mol_entries.pickle')
    #     os.system("rm %s" % network_folder + '/unfiltered_species_report.tex')
    if not os.path.exists(network_folder + '/mol_entries.pickle'):
        mol_entries = species_filter(
            database_entries,
            mol_entries_pickle_location=network_folder + '/mol_entries.pickle',
            species_report=network_folder + '/unfiltered_species_report.tex',
            species_decision_tree=species_decision_tree,
            coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy),
        )
    else:
        with open(network_folder + '/mol_entries.pickle', "rb") as pickle_file:
            mol_entries = pickle.load(pickle_file)
        logger.info("len(mol_entries):%s" % len(mol_entries))

    # add smiles info

    os.system(f"mkdir {network_folder}/xyz")
    os.system(f"mkdir {network_folder}/mol")
    if not os.path.exists(f"{network_folder}/mol_entries_smiles.csv"):
        mol_entries_smiles_df_dict = {"libe_id": [], "smiles": [], "mol_path": []}
        for molecule_data in mol_entries:
            libe_id = molecule_data.entry_id
            formula = molecule_data.formula
            molecule = molecule_data.molecule

            formula_path = formula.replace(" ", "_")
            if not os.path.exists(f"{network_folder}/mol/{formula_path}"):
                os.system(f"mkdir {network_folder}/mol/{formula_path}")
            if not os.path.exists(f"{network_folder}/xyz/{formula_path}"):
                os.system(f"mkdir {network_folder}/xyz/{formula_path}")
            xyz_path = f"{network_folder}/xyz/{formula_path}/{libe_id}.xyz"
            mol_path = f"{network_folder}/mol/{formula_path}/{libe_id}.mol"
            molecule.to(filename=xyz_path)
            molecule.to(filename=mol_path)
            mol = AllChem.MolFromMolFile(mol_path)
            try:
                Chem.SanitizeMol(mol)
                AllChem.RemoveStereochemistry(mol)
                smiles = AllChem.MolToSmiles(mol)
                # os.system(f"mv {network_folder}/{libe_id}.mol {network_folder}/{smiles}_{libe_id}.mol")
            except:
                smiles = ""
            mol_entries_smiles_df_dict["libe_id"].append(libe_id)
            mol_entries_smiles_df_dict["smiles"].append(smiles)
            mol_entries_smiles_df_dict["mol_path"].append(mol_path)

        mol_entries_smiles_df = pandas.DataFrame(mol_entries_smiles_df_dict)
        mol_entries_smiles_df.to_csv(f"{network_folder}/mol_entries_smiles.csv")

    libe_id_dict = {}
    mol_entries_smiles_df = pandas.read_csv(f"{network_folder}/mol_entries_smiles.csv")
    for index, row in mol_entries_smiles_df.iterrows():
        libe_id_dict[row["libe_id"]] = {"smiles": row["smiles"], "mol_path": row["mol_path"]}

    for molecule_data in mol_entries:
        libe_id = molecule_data.entry_id
        smiles = libe_id_dict[libe_id]["smiles"]
        molecule_data.smiles = smiles

    # Once we have generated our molecule list, we generate the bucket database
    # which is how we break up the reaction filtering amongst all avaliable workers.
    # It gets stored in the buckets.sqlite database.
    if not os.path.exists(network_folder + '/buckets.sqlite'):
        bucket(mol_entries, network_folder + '/buckets.sqlite')

    # Reaction filtering is paralellized using MPI, so we need to spawn
    # an MPI instance to run it. This is why we can't just start
    # reaction filtering by calling a python function. We pass the
    # reaction decision tree, the logging decision tree, and the electron
    # free energy as strings across this barrier. Every possible
    # reaction gets passed through both the reaction decision tree and
    # the logging decision tree. If a reaction passes the reaction
    # decision tree, it gets written to the network. If a reaction
    # passes the logging decision tree, it gets logged to the reaction
    # report along with what happened to it in reaction_decision_tree.

    # The reaction decision trees are constructed in
    # HiPRGen.reaction_questions

    params = {
        'temperature': ROOM_TEMP,
        'electron_free_energy': -1.4
    }

    dispatcher_payload = DispatcherPayload(
        network_folder + '/buckets.sqlite',
        network_folder + '/rn.sqlite',
        network_folder + '/reaction_report.tex'
    )

    worker_payload = WorkerPayload(
        network_folder + '/buckets.sqlite',
        default_reaction_decision_tree,
        params,
        Terminal.DISCARD
    )

    # The dispatcher and worker payloads are passed through the MPI barrier
    # as JSON blobs dispatcher_payload and worker_payload
    if not os.path.exists(network_folder + '/dispatcher_payload.json'):
        dumpfn(dispatcher_payload, network_folder + '/dispatcher_payload.json')
    else:
        logger.info(network_folder + '/dispatcher_payload.json' + "exists")
    if not os.path.exists(network_folder + '/worker_payload.json'):
        dumpfn(worker_payload, network_folder + '/worker_payload.json')
    else:
        logger.info(network_folder + '/worker_payload.json'+ "exists")

    if not os.path.exists(network_folder + '/rn.sqlite'):
        subprocess.run(
            [
                'mpirun',
                '--use-hwthread-cpus',
                '-n',
                number_of_threads,
                'python',
                'run_network_generation.py',
                network_folder + '/mol_entries.pickle',
                network_folder + '/dispatcher_payload.json',
                network_folder + '/worker_payload.json'
            ]
        )
    else:
        logger.info(network_folder + '/rn.sqlite' + " exists")

    return mol_entries


def li_run(network_json_path, network_folder, work_dir, init_molecule_list,
           observed_molecule_list, molecule_type, ):
    """
    :param network_json_path:
    :param work_dir: it must not exist
    :return:
    """
    mol_entries = gen_network(network_json_path,network_folder)

    # os.system("mkdir -p %s" % (work_dir))

    # After we have generated the mol_entries, we refer to molecules by
    # their index. The function find_mol_entry_from_xyz_and_charge can
    # help find the indices of specific species to be used in the initial
    # condition for propagating trajectories and/or trajectory analysis.
    Li_plus_id = find_mol_entry_from_xyz_and_charge(
        mol_entries,
        '/root/HiPRGen/xyz_files/Li.xyz',
        1)
    logger.info("Li_plus_id:%s"%Li_plus_id)
    logger.info("type(Li_plus_id):%s"%type(Li_plus_id))

    initial_state = {}

    for an_input_smiles in init_molecule_list:
        smiles_to_xyz(smiles=an_input_smiles, output_file_path='input.xyz')
        input_molecule_index = find_mol_entry_from_xyz_and_charge(mol_entries=mol_entries,
                                                                  xyz_file_path='input.xyz',
                                                                  charge=0)
        if input_molecule_index == None:
            logger.info(f'Cannot find {an_input_smiles}')
            continue
        initial_state[input_molecule_index] = 30

    if initial_state == {}:
        raise Exception("initial_state=={}")

    initial_state[Li_plus_id] = 30

    logger.info("initial_state:%s"%initial_state)
    initial_state_path = os.path.join(work_dir, 'initial_state.sqlite')
    if os.path.exists(initial_state_path):
        os.remove(initial_state_path)

    # The initial state and the trajectories (after simulation) are stored in
    # a seperate database from the network, here called initial_state.sqlite.
    # This facilitates running multiple independent simulations of the same
    # network with different initial conditions at the same time, if desired.
    insert_initial_state(initial_state, mol_entries, initial_state_path)

    # GMC is a high performance reaction network Monte Carlo simulator using the
    # Gillespie algorithm: https://github.com/BlauGroup/RNMC. Here we run 1000
    # trajectories each of 200 steps.

    logger.info("number_of_threads:%s"%number_of_threads)

    reaction_database = os.path.join(network_folder, 'rn.sqlite')
    initial_state_database = os.path.join(work_dir, 'initial_state.sqlite')

    #number_of_threads = str(16)
    command_line = f'GMC --reaction_database={reaction_database} --initial_state_database={initial_state_database} --number_of_simulations=1000 --base_seed=1000 --thread_count={number_of_threads} --step_cutoff=200'


    logger.info("command_line:%s"%command_line)

    os.system(command_line)

    logger.info('Here again')

    logger.info("WORKDIR:"% work_dir )
    os.system("ls %s" % (work_dir))
    # The network loader builds a python object around a reaction network
    # and the molecules to make it easier to use them.
    network_loader = NetworkLoader(
        network_folder + '/rn.sqlite',
        network_folder + '/mol_entries.pickle',
        work_dir + '/initial_state.sqlite'
    )

    network_loader.load_trajectories()
    network_loader.load_initial_state()

    # HiPRGen has analysis tools to understand what happened in our simulation.
    # The output files are written into the same folder in which the reaction
    # network is stored.

    # This report is empty, but we use it to generate the molecule pictures.
    # This is an expensive operation, so we only want do do it once.
    report_generator = ReportGenerator(
        network_loader.mol_entries,
        work_dir + '/dummy.tex',
        rebuild_mol_pictures=True)

    # The tally report shows reactions sorted by the number of times fired.
    reaction_tally_report(
        network_loader,
        work_dir + '/reaction_tally.tex'
    )
    os.system("cd %s && pdflatex reaction_tally.tex " % (work_dir,))

    # Run `pdflatex reaction_tally.tex` in `scratch/li_test` to generate
    # the tally report PDF.
    logger.info("WORKDIR:"% work_dir )
    os.system("ls %s" % (work_dir))
    # The species report shows every specie in the network and their IDs.
    species_report(network_loader, work_dir + '/species_report.tex')
    os.system("cd %s && pdflatex species_report.tex " % (work_dir,))
    # Run `pdflatex species_report.tex` in `scratch/li_test` to generate
    # the species report PDF.

    # Pathfinding is a central goal of HiPRGen / GMC. See mc_analysis.py for
    # further documentation of the Pathfinding class.
    pathfinding = Pathfinding(network_loader)

    # The simulation replayer sweeps through all trajectories in order
    # to extract additional information that is used for consumption
    # reports and sink reports.
    simulation_replayer = SimulationReplayer(network_loader)

    # The sink report shows species which have a production to
    # consumption ratio of greater than 3/2 and which have an expected
    # value above 0.1. These are two of the three heuristic criteria
    # that we use to identify network products. The third criteria is
    # that each network product must have a shortest path with cost
    # less than 10. This can be checked by generating pathway reports
    # to each species shown in the sink report. For the curious reader,
    # we note that generating pathway reports to the six species in the
    # sink report will show that only Li2CO3, C2H4, LiEDC-, and DLEMC
    # have sufficiently low-cost paths to pass the third criteria and
    # thus to be considered products of the test network used here.
    sink_report(simulation_replayer, work_dir + '/sink_report.tex')
    os.system("cd %s && pdflatex sink_report.tex " % (work_dir,))
    # Run `pdflatex sink_report.tex` in `scratch/li_test` to generate
    # the sink report PDF.
    logger.info("WORKDIR:"% work_dir )
    os.system("ls %s" % (work_dir))

    for an_observed_molecule in observed_molecule_list:
        smiles_to_xyz(smiles=an_observed_molecule, output_file_path='observe.xyz')
        observed_molecule_index = find_mol_entry_from_xyz_and_charge(mol_entries=mol_entries,
                                                                  xyz_file_path='observe.xyz',
                                                                  charge=0)

        if observed_molecule_index == None:
            logger.info(f'Cannot find {an_observed_molecule}')
            continue

        generate_pathway_report(
            pathfinding,
            observed_molecule_index,
            work_dir + '/observed_molecule_%s_pathways.tex' % (observed_molecule_index),
            sort_by_frequency=False)
        os.system("cd %s && pdflatex observed_molecule_%s_pathways.tex " % (work_dir, observed_molecule_index,))

        consumption_report(simulation_replayer,
                           observed_molecule_index,
                           work_dir + '/observed_molecule_%s_consumption_report.tex' % (observed_molecule_index))
        os.system("cd %s && pdflatex observed_molecule_%s_consumption_report.tex " % (work_dir, observed_molecule_index,))


if __name__ == '__main__':
    from argparse import ArgumentParser
    import json
    import rdkit
    from rdkit.Chem import AllChem
    from rdkit import Chem

    parser = ArgumentParser()

    parser.add_argument('--network_json_path',
                        type=str)  # "/root/HiPRGen/data/all_flicho_test.json"  "/root/HiPRGen/data/choli_v4.json"
    parser.add_argument('--work_dir', type=str)
    parser.add_argument('--molecule_type',
                        type=str)  # smiles or SMILES or libe_id or find_formula #EC:"O=C1OCCO1" libe-115834 LEDC: libe-115795
    parser.add_argument('--number_of_threads', default=None, type=int)
    # if molecule_type="libe_id"
    parser.add_argument('--init_molecule_libe_id_list', default="", type=str)  # EC_id libe-115834
    parser.add_argument('--observed_molecule_libe_id_list', default="", type=str)  # LEDC_id libe-115795
    # if molecule_type="SMILES" or "smiles"
    parser.add_argument('--init_molecule_smiles_list_file', default="", type=str)  #
    parser.add_argument('--observed_molecule_smiles_list_file', default="", type=str)  #

    # if molecule_type=="find_formula"
    parser.add_argument('--find_libe_id_by_formula_alphabetical_list', default="", type=str)  # split by ,

    args = parser.parse_args()

    work_dir = args.work_dir

    network_json_path = args.network_json_path
    network_json_dir, network_json_filename = os.path.split(network_json_path)
    network_folder_name = network_json_filename.split(".")[0]
    network_folder = os.path.join(network_json_dir, network_folder_name)

    number_of_threads = str(args.number_of_threads)
    if args.molecule_type not in ["SMILES", "libe_id", "smiles", "find_formula"]:
        raise Exception("""args.molecule_type not in ["SMILES","libe_id","smiles","find_formula"]""")

    if args.molecule_type == "libe_id":
        if len(args.init_molecule_libe_id_list) == 0:
            with open(f"{work_dir}/ERROR.log", "w") as fp:
                print(" len(args.init_molecule_libe_id_list)==0", file=fp)
            raise Exception(" len(args.init_molecule_libe_id_list)==0")
        init_molecule_libe_id_list = args.init_molecule_libe_id_list.split(",")

        if len(args.observed_molecule_libe_id_list) == 0:
            observed_molecule_libe_id_list = []
        else:
            observed_molecule_libe_id_list = args.observed_molecule_libe_id_list.split(",")

        li_run(network_json_path=network_json_path,network_folder=network_folder, work_dir=work_dir,
               init_molecule_list=init_molecule_libe_id_list,
               observed_molecule_list=observed_molecule_libe_id_list, molecule_type="libe_id")

    elif args.molecule_type == "SMILES" or args.molecule_type == "smiles":

        if args.init_molecule_smiles_list_file=="":
            logger.error(" len(args.init_molecule_smiles_list_file)==0")
            raise Exception("len(args.init_molecule_smiles_list)==0")
        with open(args.init_molecule_smiles_list_file, "r") as fp:

            init_molecule_smiles_list = []
            line = fp.read()
            for smiles in line.split("."):
                mol = AllChem.MolFromSmiles(smiles)
                if not mol:
                    logger.error(f"init_molecule_smiles_list :{smiles} error" )
                    raise Exception(f"init_molecule_smiles_list :{smiles} error ")
            init_molecule_smiles_list.append(AllChem.MolToSmiles(mol))


        if  args.observed_molecule_smiles_list_file == "":
            logger.info(" no observed_molecule_smiles_list_file ")
            observed_molecule_smiles_list = []
        else:
            with open(args.observed_molecule_smiles_list_file, "r") as fp:

                observed_molecule_smiles_list = []
                line = fp.read()
                for smiles in line.split("."):
                    mol = AllChem.MolFromSmiles(smiles)
                    if not mol:
                        logger.error(f"observed_molecule_smiles_list :{smiles} error")
                        raise Exception(f"observed_molecule_smiles_list :{smiles} error ")
                observed_molecule_smiles_list.append(AllChem.MolToSmiles(mol))
        logger.info(f"init_molecule_smiles_list:{init_molecule_smiles_list}")
        logger.info(f"observed_molecule_smiles_list:{observed_molecule_smiles_list}")
        li_run(network_json_path=network_json_path, network_folder=network_folder, work_dir=work_dir,
               init_molecule_list=init_molecule_smiles_list,
               observed_molecule_list=observed_molecule_smiles_list, molecule_type="smiles")
    elif args.molecule_type == "find_formula":
        formula_list = args.find_libe_id_by_formula_alphabetical_list.strip().split(",")
        with open(f"{work_dir}/EXISTING_FORMULA", "w") as fp:
            print("\n".join(os.listdir(f"{network_folder}/mol")), file=fp)
        for formula in formula_list:
            formula_path = formula.replace(" ", "_")
            formula_mol_folder_path = f"{network_folder}/mol/{formula_path}"
            formula_xyz_folder_path = f"{network_folder}/xyz/{formula_path}"
            if os.path.exists(formula_mol_folder_path):
                os.system(f"cp -r {formula_mol_folder_path} {work_dir}")
            if os.path.exists(formula_xyz_folder_path):
                os.system(f"cp -r {formula_xyz_folder_path} {work_dir}")
