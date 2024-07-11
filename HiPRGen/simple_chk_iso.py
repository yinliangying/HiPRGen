import ase
from ase.db import connect
from ase.io import read
from ase.neighborlist import build_neighbor_list
import networkx as nx
from networkx import isomorphism
import numpy as np


def get_G(atoms: ase.Atoms):
    nb_list = build_neighbor_list(atoms, self_interaction=False)
    mat = nb_list.get_connectivity_matrix()
    adj = mat.toarray()
    G = nx.from_numpy_array(adj)
    return G


def test_iso_for_two_graphs(g_1: nx.Graph, g_2: nx.Graph):
    gm = isomorphism.GraphMatcher(g_1, g_2)
    return gm.is_isomorphic()


# Read the meta.xyz file
meta_atoms = read('ortho.xyz')

# Get the charge of meta_atoms (assuming it's stored in meta_atoms.info['charge'])
# If it's stored differently, adjust this line accordingly
meta_charge = 0

if meta_charge is None:
    print("Warning: Charge information not found in meta.xyz file.")

# Calculate the graph for meta_atoms once
meta_graph = get_G(meta_atoms)

match_found = False

# Connect to the database and iterate through its entries
with connect('dump.db') as db:
    for row in db.select():
        numbers_identical = np.array_equal(row.numbers, meta_atoms.numbers)
        if not numbers_identical:
            continue
        db_atoms = row.toatoms()
        db_charge = row.data.get('charge', None)

        # Check if charges are equal (if charge information is available)
        charges_equal = (meta_charge == db_charge) if meta_charge is not None and db_charge is not None else None

        # Calculate graph for db_atoms and check if topologies are isomorphic
        db_graph = get_G(db_atoms)
        topo_isomorphic = test_iso_for_two_graphs(meta_graph, db_graph)

        if topo_isomorphic and charges_equal:
            print(f"Match found! Row id: {row.id}")
            print(f"Topology: Isomorphic, Charge: {db_charge}")
            match_found = True
            break  # Exit the loop after finding the first match

if not match_found:
    print("No matching atoms found in the database.")