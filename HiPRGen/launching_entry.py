from pdf2image import convert_from_path
import os
import pickle
import sqlite3

from HiPRGen.initial_state import insert_initial_state, find_mol_entry_from_xyz_and_charge
from HiPRGen.mc_analysis import (
    reaction_tally_report,
    Pathfinding,
    SimulationReplayer,
    generate_pathway_report,
    sink_report,
    consumption_report
)
from HiPRGen.network_loader import NetworkLoader
from rdkit import Chem
from rdkit.Chem import AllChem


# libe_default_paths = {
#     'rn_db_path': r'/root/HiPRGen/data/libe/rn.sqlite',
#     'mol_entry_file_path': r'/root/HiPRGen/data/libe/mol_entries.pickle',
#     'mol_picture_folder_path': r'/root/HiPRGen/data/libe/mol_pictures',
# }
database_dir= "/root/HiPRGen/data/new_libe_fmol_20240731/"
libe_default_paths = {
    'rn_db_path':rf'{database_dir}/rn.sqlite',
    'mol_entry_file_path':rf'{database_dir}/mol_entries.pickle',
    'mol_picture_folder_path':rf'{database_dir}/mol_pictures',
}


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


def tex_to_results(tex_file_path):
    cwd_tex = os.getcwd()
    abs_tex_file_path = os.path.abspath(tex_file_path)
    workbase, tex_file_name = os.path.split(abs_tex_file_path)
    os.chdir(workbase)
    os.system(f'pdflatex {tex_file_name}')
    for a_file in os.listdir('./'):
        if a_file.endswith('.tex') or a_file.endswith('.aux') or a_file.endswith('.log'):
            os.remove(a_file)
    pure_file_name = os.path.splitext(tex_file_name)[0]
    images = convert_from_path(f'{pure_file_name}.pdf')
    png_folder_name = pure_file_name + '_png'
    os.makedirs(png_folder_name)
    os.chdir(png_folder_name)
    for i, image in enumerate(images):
        image.save(f"page_{i + 1}.png", "PNG")
    os.chdir(cwd_tex)


def chk_ids_with_rn_db(main_mol_id: int, sub_mol_id: int, rn_db_path: str):
    conn = sqlite3.connect(rn_db_path)
    cursor = conn.cursor()
    query = """
      SELECT EXISTS (
        SELECT 1 FROM reactions
        WHERE (reactant_1 = ? AND reactant_2 = ?)
        OR (reactant_1 = ? AND reactant_2 = ?)
      );
    """
    cursor.execute(query, (main_mol_id, sub_mol_id, sub_mol_id, main_mol_id))
    has_reaction = cursor.fetchone()[0]
    conn.close()
    return has_reaction > 0


def run_with_id(main_mol_id: int, sub_mol_id: int, output_dir: str, simulation_times: int, num_cores: int,
                default_file_paths: dict = libe_default_paths):
    cwd_ = os.getcwd()
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    mol_entry_file_path = default_file_paths['mol_entry_file_path']
    mol_picture_folder_path = default_file_paths['mol_picture_folder_path']
    rn_db_path = default_file_paths['rn_db_path']
    match_flag = chk_ids_with_rn_db(main_mol_id=main_mol_id, sub_mol_id=sub_mol_id, rn_db_path=rn_db_path)
    if not match_flag:
        os.chdir(output_dir)
        with open('No match reaction id for input', 'w') as f:
            pass
        os.chdir(cwd_)
        return

    with open(mol_entry_file_path, "rb") as pickle_file:
        mol_entries = pickle.load(pickle_file)
    initial_state = {main_mol_id: 30, sub_mol_id: 30}
    initial_state_path = os.path.join(output_dir, 'initial_state.sqlite')
    insert_initial_state(initial_state, mol_entries, initial_state_path)

    command_line = f'GMC --reaction_database={rn_db_path} --initial_state_database={initial_state_path} --number_of_simulations={simulation_times} --base_seed=1000 --thread_count={num_cores} --step_cutoff=200'
    os.system(command_line)

    network_loader = NetworkLoader(
        network_database=rn_db_path,
        mol_entries_pickle=mol_entry_file_path,
        initial_state_database=initial_state_path
    )
    network_loader.load_trajectories()
    network_loader.load_initial_state()
    os.chdir(output_dir)
    reaction_tally_report(network_loader=network_loader, reaction_tally_report_path='reaction_tally.tex',
                          fixed_mol_pictures_folder=mol_picture_folder_path)
    tex_to_results(tex_file_path='reaction_tally.tex')
    simulation_replayer = SimulationReplayer(network_loader=network_loader)
    sink_report(simulation_replayer=simulation_replayer, sink_report_path='sink_report.tex',
                fixed_mol_pictures_folder=mol_picture_folder_path)
    tex_to_results(tex_file_path='sink_report.tex')
    # get sink mol id and dump pathways
    pathfinding = Pathfinding(network_loader)
    simulation_replayer = SimulationReplayer(network_loader=network_loader)
    sinks_sorted = sorted(
        simulation_replayer.sinks,
        key=lambda i: -simulation_replayer.sink_data[i]["ratio"])
    pathways_abs_workbase = os.path.abspath('pathways')
    consumption_abs_workbase = os.path.abspath('consumption')
    os.makedirs(pathways_abs_workbase)
    os.makedirs(consumption_abs_workbase)
    for species_index in sinks_sorted[:10]:
        mol_entry = mol_entries[species_index]
        formula = mol_entry.molecule.composition.alphabetical_formula
        formula = '_'.join(formula.split())
        tex_file_name = f'id_{species_index}_formula_{formula}.tex'
        os.chdir(pathways_abs_workbase)
        generate_pathway_report(
            pathfinding=pathfinding,
            species_id=species_index,
            report_file_path=tex_file_name,
            number_of_pathways=10,
            fixed_mol_pictures_folder=mol_picture_folder_path
        )
        tex_to_results(tex_file_path=tex_file_name)
        os.chdir(consumption_abs_workbase)
        consumption_report(
            simulation_replayer=simulation_replayer,
            species_index=species_index,
            consumption_report_path=tex_file_name,
            fixed_mol_pictures_folder=mol_picture_folder_path
        )
        tex_to_results(tex_file_path=tex_file_name)
    os.chdir(cwd_)


def run_with_smiles(output_dir: str, simulation_times: int, num_cores: int, main_mol_smiles: str, sub_mol_name: str=None,
                default_file_paths: dict = libe_default_paths):
    mol_entry_file_path = default_file_paths['mol_entry_file_path']
    with open(mol_entry_file_path, "rb") as pickle_file:
        mol_entries = pickle.load(pickle_file)

    cwd_ = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)
    try:
        smiles_to_xyz(smiles=main_mol_smiles, output_file_path=r'main_mol.xyz')
        main_mol_id = find_mol_entry_from_xyz_and_charge(mol_entries=mol_entries,
                                                         xyz_file_path=r'main_mol.xyz',
                                                         charge=0)
    except:
        with open('No species match with the input SMILES', 'w') as f:
            pass
        return
    os.chdir(cwd_)
    if sub_mol_name == None:
        sub_id = -1
    else:
        raise NotImplementedError
        common_sub_mol_info = {'OH-1': 1835, 'C2H4': 9998, 'H2O': 7689, 'H2': 8930, 'CO': 4293, 'CO2': 579, 'H+1': 8450,
                               'Li2CO3+1': 8380, 'LiCO3-1': 7950, 'Li+1': 5253, 'LiF': 5259} #这个是旧libe的id
        sub_id = common_sub_mol_info[sub_mol_name]
    run_with_id(main_mol_id=main_mol_id, sub_mol_id=sub_id, output_dir=output_dir,
                simulation_times=simulation_times, num_cores=num_cores, default_file_paths=default_file_paths)
