import math
from HiPRGen.mol_entry import MoleculeEntry
from functools import partial
import itertools
import networkx as nx
from networkx.algorithms.graph_hashing import weisfeiler_lehman_graph_hash
from HiPRGen.constants import Terminal, ROOM_TEMP, KB, PLANCK, m_formulas
from monty.json import MSONable

"""
The reaction decision tree:

A question is a function q(reaction, mol_entries, params) -> Bool

reaction is a dict:

        reaction = {
            'reactants' : reactant indices
            'products' : product indices,
            'number_of_reactants',
            'number_of_products'}

params is a dict:


        params = {
           'temperature'
           'electron_free_energy'
        }

The lists of reactant and product indices always have length two. We use -1 when there is a only a single reactant or product.

The questions can also set reaction['rate'] and reaction['dG']

Questions will be writable by hand, or we could have machine learning filters.

A node is either a Terminal or a non empty list [(question, node)]

class Terminal(Enum):
    KEEP = 1
    DISCARD = -1

For the return value of a question, True means travel to this node and False means try next question in the list.

for non terminal nodes, it is an error if every question returns False. i.e getting stuck at a non terminal node is an error.

Once a Terminal node is reached, it tells us whether to keep or discard the reaction.

logging decision tree: The dispatcher takes a second decision tree as an argument, the logging decision tree. Reactions which return Terminal.KEEP from the logging decision tree will be logged in the generation report, with location specified by the argument generation_report_path

"""

hydrogen_graph = nx.MultiGraph()
hydrogen_graph.add_node(0, specie='H')
hydrogen_hash = weisfeiler_lehman_graph_hash(
    hydrogen_graph,
    node_attr='specie')


def run_decision_tree(
        reaction,
        mol_entries,
        params,
        decision_tree,
        decision_pathway=None):
    node = decision_tree

    while type(node) == list:
        next_node = None
        for (question, new_node) in node:
            if question(reaction, mol_entries, params):

                # if decision_pathway is a list,
                # append the question which
                # answered true i.e the edge we follow
                if decision_pathway is not None:
                    decision_pathway.append(question)

                next_node = new_node
                break

        node = next_node


    if type(node) == Terminal:
        if decision_pathway is not None:
            decision_pathway.append(node)

        if node == Terminal.KEEP:
            return True
        else:
            return False
    else:
        print(node)
        raise Exception(
            """
            unexpected node type reached.
            this is usually caused because none of the questions in some node returned True.
            """)



def default_rate(dG_barrier, params):
    kT = KB * params['temperature']
    max_rate = kT / PLANCK
    rate = max_rate * math.exp(- dG_barrier / kT)
    return rate

class dG_above_threshold(MSONable):

    def __init__(self, threshold, free_energy_type, constant_barrier):

        self.threshold = threshold
        self.free_energy_type = free_energy_type
        self.constant_barrier = constant_barrier

        if free_energy_type == 'free_energy':
            self.get_free_energy = lambda mol: mol.free_energy
        elif free_energy_type == 'solvation_free_energy':
            self.get_free_energy = lambda mol: mol.solvation_free_energy
        else:
            raise Exception("unrecognized free energy type")


    def __call__(self, reaction, mol_entries, params):


        dG = 0.0

        # positive dCharge means electrons are lost
        dCharge = 0.0

        for i in range(reaction['number_of_reactants']):
            reactant_index = reaction['reactants'][i]
            mol = mol_entries[reactant_index]
            dG -= self.get_free_energy(mol)
            dCharge -= mol.charge

        for j in range(reaction['number_of_products']):
            product_index = reaction['products'][j]
            mol = mol_entries[product_index]
            dG += self.get_free_energy(mol)
            dCharge += mol.charge

        dG += dCharge * params['electron_free_energy']

        if dG > self.threshold:
            return True
        else:
            reaction['dG'] = dG
            if dG < 0:
                barrier = self.constant_barrier
            else:
                barrier = dG + self.constant_barrier

            reaction['dG_barrier'] = barrier
            reaction['rate'] = default_rate(barrier, params)
            return False



class is_redox_reaction(MSONable):

    def __init__(self):
        pass

    def __call__(self, reaction, mol_entries, params):
        # positive dCharge means electrons are lost
        dCharge = 0.0

        for i in range(reaction['number_of_reactants']):
            reactant_index = reaction['reactants'][i]
            mol = mol_entries[reactant_index]
            dCharge -= mol.charge

        for j in range(reaction['number_of_products']):
            product_index = reaction['products'][j]
            mol = mol_entries[product_index]
            dCharge += mol.charge

        if dCharge == 0:
            return False
        else:
            return True


class too_many_reactants_or_products(MSONable):
    def __init__(self):
        pass


    def __call__(self, reaction, mols, params):
        if (reaction['number_of_reactants'] != 1 or
            reaction['number_of_products'] != 1):
            return True
        else:
            return False




class dcharge_too_large(MSONable):
    def __init__(self):
        pass

    def __call__(self, reaction, mol_entries, params):
        dCharge = 0.0

        for i in range(reaction['number_of_reactants']):
            reactant_index = reaction['reactants'][i]
            mol = mol_entries[reactant_index]
            dCharge -= mol.charge

        for j in range(reaction['number_of_products']):
            product_index = reaction['products'][j]
            mol = mol_entries[product_index]
            dCharge += mol.charge

        if abs(dCharge) > 1:
            return True
        else:
            return False



class set_redox_rate_markus_theory(MSONable):
    """
        Okay, so Marcus Theory.The math works out like so.∆G* = λ/4 (1 +
    ∆G / λ)^2 ∆G is the Gibbs free energy of the reaction, ∆G* is the
    energy barrier, and λ is the “reorganization energy” (basically the
    energy penalty for reorganizing the solvent environment to accommodate
    the change in local charge).The reorganization energy can be broken up
    into two terms, an inner term (“i”) representing the contribution from
    the first solvation shell and an outer term (“o”) representing the
    contribution from the bulk solvent: λ = λi + λoλo = ∆e/(8 pi ε0) (1/r
    - 1/R) (1/n^2 - 1/ε) where ∆e is the change in charge in terms of
    fundamental charge (1.602 * 10 ^-19 C), ε0 is the vacuum permittivity
    (8.854 * 10 ^-12 F/m), r is the first solvation shell radius (I
    usually just pick a constant, say 6 Angstrom), R is the distance to
    the electrode (again, for these purposes, just pick something - say
    7.5 Angstrom), n is the index of refraction (1.415 for EC) and ε is
    the relative dielectric (18.5 for EC/EMC).
    """

    def __init__(self, threshold):
        self.threshold = threshold
        pass



    def __call__(self, reaction, mol_entries, params):
        dCharge = 0.0

        reactant_index = reaction['reactants'][0]
        reactant = mol_entries[reactant_index]
        dCharge -= reactant.charge

        product_index = reaction['products'][0]
        product = mol_entries[product_index]
        dCharge += product.charge

        dG_barrier = float("inf")

        try:
            product_coordimers = product.coordimers.values()
        except:
            product_coordimers = [product]



        for m in product_coordimers:

            n = 1.415  # index of refraction; variable
            eps = 18.5  # dielectric constant; variable

            r = 6.0  # in Angstrom
            R = 7.5  # in Angstrom

            eps_0 = 8.85419 * 10 ** -12  # vacuum permittivity
            e = 1.602 * 10 ** -19  # fundamental charge

            l_outer = e / (8 * math.pi * eps_0)
            l_outer *= (1 / r - 1/(2 * R)) * 10 ** 10  # Converting to SI units; factor of 2 is because of different definitions of the distance to electrode
            l_outer *= (1 / n ** 2 - 1 / eps)

            if dCharge == -1:
                vals = [reactant.electron_affinity, m.ionization_energy]
                vals_filtered = [v for v in vals if v is not None]
                l_inner = sum(vals_filtered) / len(vals_filtered)

            if dCharge == 1:
                vals = [reactant.ionization_energy, m.electron_affinity]
                vals_filtered = [v for v in vals if v is not None]
                l_inner = sum(vals_filtered) / len(vals_filtered)


            if l_inner < 0:
                l_inner = 0

            l = l_inner + l_outer


            dG_temp = m.free_energy - reactant.free_energy + dCharge * params['electron_free_energy']
            dG_barrier_temp = l / 4 * (1 + dG_temp / l) ** 2

            dG_barrier = min(dG_barrier, dG_barrier_temp)


        if dG_barrier < 0:
            dG_barrier = 0

        dG = product.free_energy - reactant.free_energy + dCharge * params['electron_free_energy']

        reaction['dG'] = dG
        reaction['dG_barrier'] = dG_barrier
        reaction['rate'] = default_rate(dG_barrier, params)

        if dG > self.threshold:
            return True
        else:
            return False


class reactant_and_product_not_isomorphic(MSONable):

    def __init__(self):
        pass

    def __call__(self, reaction, mols, params):
        reactant = mols[reaction['reactants'][0]]
        product = mols[reaction['products'][0]]
        if reactant.covalent_hash != product.covalent_hash:
            return True
        else:
            return False


class default_true(MSONable):

    def __init__(self):
        pass

    def __call__(self, reaction, mols, params):
        return True

def star_count_diff_above_threshold(threshold, reaction, mols, params):
    reactant_stars = {}
    product_stars = {}
    tags = set()

    for i in range(reaction['number_of_reactants']):
        reactant_index = reaction['reactants'][i]
        mol = mols[reactant_index]
        for h in mol.star_hashes.values():
            tags.add(h)
            if h in reactant_stars:
                reactant_stars[h] += 1
            else:
                reactant_stars[h] = 1

    for j in range(reaction['number_of_products']):
        product_index = reaction['products'][j]
        mol = mols[product_index]
        for h in mol.star_hashes.values():
            tags.add(h)
            if h in product_stars:
                product_stars[h] += 1
            else:
                product_stars[h] = 1

    count = 0

    for tag in tags:
        count += abs(reactant_stars.get(tag,0) - product_stars.get(tag,0))

    if count > threshold:
        return True
    else:
        return False

def reaction_is_covalent_decomposable(reaction, mols, params):
    if (reaction['number_of_reactants'] == 2 and
        reaction['number_of_products'] == 2):


        reactant_total_hashes = set()
        for i in range(reaction['number_of_reactants']):
            reactant_id = reaction['reactants'][i]
            reactant = mols[reactant_id]
            reactant_total_hashes.add(reactant.covalent_hash)

        product_total_hashes = set()
        for i in range(reaction['number_of_products']):
            product_id = reaction['products'][i]
            product = mols[product_id]
            product_total_hashes.add(product.covalent_hash)

        if len(reactant_total_hashes.intersection(product_total_hashes)) > 0:
            return True
        else:
            return False

    return False


def metal_coordination_passthrough(reaction, mols, params):

    for i in range(reaction['number_of_reactants']):
        reactant_id = reaction['reactants'][i]
        reactant = mols[reactant_id]
        if reactant.formula in m_formulas:
            return True

    for i in range(reaction['number_of_products']):
        product_id = reaction['products'][i]
        product = mols[product_id]
        if product.formula in m_formulas:
            return True

    return False



def fragment_matching_found(reaction, mols, params):

    reactant_fragment_indices_list = []
    product_fragment_indices_list = []

    if reaction['number_of_reactants'] == 1:
        reactant = mols[reaction['reactants'][0]]
        for i in range(len(reactant.fragment_data)):
            reactant_fragment_indices_list.append([i])


    if reaction['number_of_reactants'] == 2:
        reactant_0 = mols[reaction['reactants'][0]]
        reactant_1 = mols[reaction['reactants'][1]]
        for i in range(len(reactant_0.fragment_data)):
            for j in range(len(reactant_1.fragment_data)):
                if (reactant_0.fragment_data[i].number_of_bonds_broken +
                    reactant_1.fragment_data[j].number_of_bonds_broken <= 1):

                    reactant_fragment_indices_list.append([i,j])


    if reaction['number_of_products'] == 1:
        product = mols[reaction['products'][0]]
        for i in range(len(product.fragment_data)):
            product_fragment_indices_list.append([i])


    if reaction['number_of_products'] == 2:
        product_0 = mols[reaction['products'][0]]
        product_1 = mols[reaction['products'][1]]
        for i in range(len(product_0.fragment_data)):
            for j in range(len(product_1.fragment_data)):
                if (product_0.fragment_data[i].number_of_bonds_broken +
                    product_1.fragment_data[j].number_of_bonds_broken <= 1):

                    product_fragment_indices_list.append([i,j])


    for reactant_fragment_indices in reactant_fragment_indices_list:
        for product_fragment_indices in product_fragment_indices_list:
            reactant_fragment_count = 0
            product_fragment_count = 0
            reactant_bonds_broken = []
            product_bonds_broken = []

            reactant_hashes = dict()
            for reactant_index, frag_complex_index in enumerate(
                    reactant_fragment_indices):

                fragment_complex = mols[
                    reaction['reactants'][reactant_index]].fragment_data[
                        frag_complex_index]

                for bond in fragment_complex.bonds_broken:
                    reactant_bonds_broken.append([(reactant_index, x) for x in bond])

                for i in range(fragment_complex.number_of_fragments):
                    reactant_fragment_count += 1
                    tag = fragment_complex.fragment_hashes[i]
                    if tag in reactant_hashes:
                        reactant_hashes[tag] += 1
                    else:
                        reactant_hashes[tag] = 1

            product_hashes = dict()
            for product_index, frag_complex_index in enumerate(
                    product_fragment_indices):

                fragment_complex = mols[
                    reaction['products'][product_index]].fragment_data[
                        frag_complex_index]

                for bond in fragment_complex.bonds_broken:
                    product_bonds_broken.append([(product_index, x) for x in bond])


                for i in range(fragment_complex.number_of_fragments):
                    product_fragment_count += 1
                    tag = fragment_complex.fragment_hashes[i]
                    if tag in product_hashes:
                        product_hashes[tag] += 1
                    else:
                        product_hashes[tag] = 1


            # don't consider fragmentations with both a ring opening and closing
            if (reaction['number_of_reactants'] == 2 and
                reaction['number_of_products'] == 2 and
                reactant_fragment_count == 2 and
                product_fragment_count == 2):
                continue


            if reactant_hashes == product_hashes:
                reaction['reactant_bonds_broken'] = reactant_bonds_broken
                reaction['product_bonds_broken'] = product_bonds_broken
                reaction['hashes'] = reactant_hashes
                reaction['reactant_fragment_count'] = reactant_fragment_count
                reaction['product_fragment_count'] = product_fragment_count

                return True

    return False

def single_reactant_single_product_not_hydrogen_transfer(reaction, mols, params):
    if (reaction['number_of_reactants'] == 1 and
        reaction['number_of_products'] == 1 and
        len(reaction['reactant_bonds_broken']) == 1 and
        len(reaction['product_bonds_broken']) == 1 and
        hydrogen_hash not in reaction['hashes']):

        return True

    return False

def single_reactant_double_product_ring_close(reaction, mols, params):
    if (reaction['number_of_reactants'] == 1 and
        reaction['number_of_products'] == 2 and
        len(reaction['reactant_bonds_broken']) == 1 and
        len(reaction['product_bonds_broken']) == 1 and
        reaction['product_fragment_count'] == 2):

        return True

    return False



def concerted_metal_coordination(reaction, mols, params):
    if (reaction['number_of_reactants'] == 2 and
        reaction['number_of_products'] == 2):

        reactant_0 = mols[reaction['reactants'][0]]
        reactant_1 = mols[reaction['reactants'][1]]
        product_0 = mols[reaction['products'][0]]
        product_1 = mols[reaction['products'][1]]



        if (reactant_0.formula in m_formulas or
            reactant_1.formula in m_formulas or
            product_0.formula in m_formulas or
            product_1.formula in m_formulas):
            return True
        else:
            return False

    return False


def concerted_metal_coordination_one_product(reaction, mols, params):

    if (reaction['number_of_reactants'] == 2 and
        reaction['number_of_products'] == 1):

        reactant_0 = mols[reaction['reactants'][0]]
        reactant_1 = mols[reaction['reactants'][1]]
        product = mols[reaction['products'][0]]

        reactant_covalent_hashes = set([
            reactant_0.covalent_hash,
            reactant_1.covalent_hash])

        if ((reactant_0.formula in m_formulas or
            reactant_1.formula in m_formulas) and
            product.covalent_hash not in reactant_covalent_hashes
            ):
            return True
        else:
            return False

    return False


def concerted_metal_coordination_one_reactant(reaction, mols, params):

    if (reaction['number_of_reactants'] == 1 and
        reaction['number_of_products'] == 2):

        product_0 = mols[reaction['products'][0]]
        product_1 = mols[reaction['products'][1]]
        reactant = mols[reaction['reactants'][0]]

        product_covalent_hashes = set([
            product_0.covalent_hash,
            product_1.covalent_hash])

        if ((product_0.formula in m_formulas or
            product_1.formula in m_formulas) and
            reactant.covalent_hash not in product_covalent_hashes
            ):
            return True
        else:
            return False

    return False



li_ec_reaction_markus_theory_decision_tree = [

    # redox branch
    (is_redox_reaction(), [

        (too_many_reactants_or_products(), Terminal.DISCARD),
        (dcharge_too_large(), Terminal.DISCARD),
        (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
        (set_redox_rate_markus_theory(0.5), Terminal.DISCARD),
        (default_true(), Terminal.KEEP)
    ]),


    (dG_above_threshold(0.5, "solvation_free_energy", 0.2), Terminal.DISCARD),

    (partial(star_count_diff_above_threshold, 4), Terminal.DISCARD),

    (reaction_is_covalent_decomposable, Terminal.DISCARD),

    (concerted_metal_coordination, Terminal.DISCARD),

    (concerted_metal_coordination_one_product, Terminal.DISCARD),

    (concerted_metal_coordination_one_reactant, Terminal.DISCARD),

    (metal_coordination_passthrough, Terminal.KEEP),

    (fragment_matching_found, [
        (single_reactant_single_product_not_hydrogen_transfer, Terminal.DISCARD),
        (single_reactant_double_product_ring_close, Terminal.DISCARD),
        (default_true(), Terminal.KEEP)]
    ),

    (default_true(), Terminal.DISCARD)
    ]


li_ec_reaction_decision_tree = [

    # redox branch
    (is_redox_reaction(), [

        (too_many_reactants_or_products(), Terminal.DISCARD),
        (dcharge_too_large(), Terminal.DISCARD),
        (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
        (dG_above_threshold(0.5, "free_energy", 0.0), Terminal.DISCARD),
        (default_true(), Terminal.KEEP)
    ]),


    (dG_above_threshold(0.5, "solvation_free_energy", 0.0), Terminal.DISCARD),

    (partial(star_count_diff_above_threshold, 4), Terminal.DISCARD),

    (reaction_is_covalent_decomposable, Terminal.DISCARD),

    (concerted_metal_coordination, Terminal.DISCARD),

    (concerted_metal_coordination_one_product, Terminal.DISCARD),

    (concerted_metal_coordination_one_reactant, Terminal.DISCARD),

    (metal_coordination_passthrough, Terminal.KEEP),

    (fragment_matching_found, [
        (single_reactant_single_product_not_hydrogen_transfer, Terminal.DISCARD),
        (single_reactant_double_product_ring_close, Terminal.DISCARD),
        (default_true(), Terminal.KEEP)]
    ),

    (default_true(), Terminal.DISCARD)
    ]



mg_g2_reaction_decision_tree = [

    # redox branch
    (is_redox_reaction(), [

        (dG_above_threshold(0.5, "free_energy", 0), Terminal.DISCARD),

        (too_many_reactants_or_products(), Terminal.DISCARD),
        (dcharge_too_large(), Terminal.DISCARD),
        (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
        (default_true(), Terminal.KEEP)
    ]),

    (dG_above_threshold(
             0.5, "solvation_free_energy", 0), Terminal.DISCARD),

    (partial(star_count_diff_above_threshold, 4), Terminal.DISCARD),

    (reaction_is_covalent_decomposable, Terminal.DISCARD),

    (concerted_metal_coordination, Terminal.DISCARD),

    (concerted_metal_coordination_one_product, Terminal.DISCARD),

    (concerted_metal_coordination_one_reactant, Terminal.DISCARD),

    (metal_coordination_passthrough, Terminal.KEEP),

    (fragment_matching_found, [
        (single_reactant_single_product_not_hydrogen_transfer, Terminal.DISCARD),
        (default_true(), Terminal.KEEP)]
    ),

    (default_true(), Terminal.DISCARD)
    ]

mg_thf_reaction_decision_tree = [

    # redox branch
    (is_redox_reaction(), [

        (dG_above_threshold(
                 0.5,
                 "free_energy", 0), Terminal.DISCARD),

        (too_many_reactants_or_products(), Terminal.DISCARD),
        (dcharge_too_large(), Terminal.DISCARD),
        (reactant_and_product_not_isomorphic(), Terminal.DISCARD),
        (default_true(), Terminal.KEEP)
    ]),

    (dG_above_threshold(
             0.5,
             "solvation_free_energy", 0), Terminal.DISCARD),

    (partial(star_count_diff_above_threshold, 4), Terminal.DISCARD),

    (reaction_is_covalent_decomposable, Terminal.DISCARD),

    (concerted_metal_coordination, Terminal.DISCARD),

    (concerted_metal_coordination_one_product, Terminal.DISCARD),

    (concerted_metal_coordination_one_reactant, Terminal.DISCARD),

    (metal_coordination_passthrough, Terminal.KEEP),

    (fragment_matching_found, [
        (single_reactant_single_product_not_hydrogen_transfer, Terminal.DISCARD),
        (default_true(), Terminal.KEEP)]
    ),

    (default_true(), Terminal.DISCARD)
    ]


standard_logging_decision_tree = Terminal.DISCARD

li_ec_redox_logging_decision_tree = Terminal.DISCARD



# this dictionary exists so that we can pass a decision tree argument to mpiexec
reaction_decision_tree_dict = {
    'li_ec_reaction_decision_tree' : li_ec_reaction_decision_tree,
    'li_ec_reaction_markus_theory_decision_tree' : li_ec_reaction_markus_theory_decision_tree,
    'mg_g2_reaction_decision_tree' : mg_g2_reaction_decision_tree,
    'mg_thf_reaction_decision_tree' : mg_thf_reaction_decision_tree
}

logging_decision_tree_dict = {
    'standard_logging_decision_tree' : standard_logging_decision_tree,
    'li_ec_redox_logging_decision_tree' : li_ec_redox_logging_decision_tree
    }

params_dict = {
    "-1.4" : {
        'temperature' : ROOM_TEMP,
        'electron_free_energy' : -1.4
    },

    "-1.9" : {
        'temperature' : ROOM_TEMP,
        'electron_free_energy' : -1.9
    },

    "-2.06" : {
        'temperature' : ROOM_TEMP,
        'electron_free_energy' : -2.06
    }
}

