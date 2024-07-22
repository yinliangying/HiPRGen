
class FragmentComplex():
    def __init__(self,number_of_bonds_broken:int=None,
                 bonds_broken:list[list[int,int]]=None,
                 number_of_fragments:int=None,
                 fragment_hashes:list[str]=None):
        self.number_of_bonds_broken=number_of_bonds_broken
        self.bonds_broken=bonds_broken
        self.number_of_fragments=number_of_fragments
        self.fragment_hashes=fragment_hashes

class MoleculeEntry():
    def __init__(self,number_of_fragment_data=None,fragment_data:list[FragmentComplex]=None):
        self.number_of_fragment_data=number_of_fragment_data
        self.fragment_data=fragment_data

class Return():
    def __init__(self,
                 r:bool=None,
                 reactant_fragment_count:int=None,product_fragment_count:int=None,
                 num_reactant_bonds_broken:int=None,num_product_bonds_broken:int=None,
                 reactant_bonds_broken:list[(int,int),(int,int)]=None,
                 product_bonds_broken:list[(int,int),(int,int)]=None,
                 hashes_key:list[str]=None,hashes_value:list[int]=None,num_hashes:int=None):
        self.r=r
        self.reactant_fragment_count=reactant_fragment_count
        self.product_fragment_count=product_fragment_count
        self.num_reactant_bonds_broken=num_reactant_bonds_broken
        self.num_product_bonds_broken=num_product_bonds_broken
        self.reactant_bonds_broken=reactant_bonds_broken
        self.product_bonds_broken=product_bonds_broken
        self.hashes_key=hashes_key
        self.hashes_value=hashes_value
        self.num_hashes=num_hashes




def function( number_of_reactants:int,number_of_products:int,
              reactant0_mol:MoleculeEntry,reactant1_mol:MoleculeEntry,
              product0_mol:MoleculeEntry,product1_mol:MoleculeEntry):

    reactant_fragment_indices_list = []
    product_fragment_indices_list = []
    result=Return()
    if number_of_reactants == 1:
        for i in range(len(reactant0_mol.fragment_data)):
            reactant_fragment_indices_list.append([i,-1])

    if number_of_reactants == 2:
        for i in range(len(reactant0_mol.fragment_data)):
            for j in range(len(reactant1_mol.fragment_data)):
                if (reactant0_mol.fragment_data[i].number_of_bonds_broken +
                        reactant1_mol.fragment_data[j].number_of_bonds_broken <= 1):
                    reactant_fragment_indices_list.append([i, j])

    if number_of_products== 1:
        for i in range(len(product0_mol.fragment_data)):
            product_fragment_indices_list.append([i,-1])

    if number_of_products == 2:
        for i in range(len(product0_mol.fragment_data)):
            for j in range(len(product1_mol.fragment_data)):
                if (product0_mol.fragment_data[i].number_of_bonds_broken +
                        product1_mol.fragment_data[j].number_of_bonds_broken <= 1):
                    product_fragment_indices_list.append([i, j])

    for reactant_fragment_indices in reactant_fragment_indices_list:
        for product_fragment_indices in product_fragment_indices_list:
            reactant_fragment_count = 0
            product_fragment_count = 0
            reactant_bonds_broken = []
            product_bonds_broken = []

            reactant_hashes = dict()
            for reactant_index, frag_complex_index in enumerate(
                    reactant_fragment_indices):

                if reactant_index == 0:
                    fragment_complex = reactant0_mol.fragment_data[
                        frag_complex_index]
                elif reactant_index == 1:
                    if frag_complex_index == -1:
                        continue
                    fragment_complex = reactant1_mol.fragment_data[
                        frag_complex_index]
                else:
                    continue


                for bond in fragment_complex.bonds_broken:
                    reactant_bonds_broken.append(
                        ((reactant_index, bond[0]), (reactant_index, bond[1]))
                    )

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


                if product_index == 0:
                    fragment_complex = product0_mol.fragment_data[
                        frag_complex_index]
                elif product_index == 1:
                    fragment_complex = product1_mol.fragment_data[
                        frag_complex_index]
                else:
                    continue

                for bond in fragment_complex.bonds_broken:
                    product_bonds_broken.append(
                        ((product_index, bond[0]), (product_index, bond[1])))

                for i in range(fragment_complex.number_of_fragments):
                    product_fragment_count += 1
                    tag = fragment_complex.fragment_hashes[i]
                    if tag in product_hashes:
                        product_hashes[tag] += 1
                    else:
                        product_hashes[tag] = 1

            # don't consider fragmentations with both a ring opening and closing
            if ( number_of_reactants == 2 and
                   number_of_products  == 2 and
                    reactant_fragment_count == 2 and
                    product_fragment_count == 2):
                continue

            if reactant_hashes == product_hashes:
                result.r=True
                result.num_hashes = len(reactant_hashes)
                result.hashes_key = list(reactant_hashes.keys())
                result.hashes_value = list(reactant_hashes.values())
                result.num_reactant_bonds_broken = len(reactant_bonds_broken)
                result.num_product_bonds_broken = len(product_bonds_broken)
                result.reactant_bonds_broken = reactant_bonds_broken
                result.product_bonds_broken = product_bonds_broken
                result.reactant_fragment_count=reactant_fragment_count
                result.product_fragment_count=product_fragment_count

                return result
    result.r=False
    return result
