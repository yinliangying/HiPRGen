
import ctypes
import pickle
from ctypes import Structure, c_int, POINTER,c_char,c_char_p,c_bool
from HiPRGen.mol_entry import MoleculeEntry
from time import time
import copy
def create_molecule_entry(mol_entries,reactant_id):
    molecule_entry_ctype = MoleculeEntry_c_type()
    if reactant_id == -1:
        return molecule_entry_ctype

    mol_entry = mol_entries[reactant_id]

    max_list_size = len(mol_entry.fragment_data)
    for f_idx, fragment_complex in enumerate(mol_entry.fragment_data):
        number_of_bonds_broken = fragment_complex.number_of_bonds_broken
        number_of_fragments = fragment_complex.number_of_fragments
        max_list_size = max(max_list_size, number_of_bonds_broken, number_of_fragments)
    if max_list_size>MAX_LIST_SIZE:
        print(f"mol_id:{reactant_id},max_list_size:{max_list_size} skip")
        return None

    number_of_fragment_data=0
    for f_idx,fragment_complex in enumerate(mol_entry.fragment_data):
        number_of_bonds_broken = fragment_complex.number_of_bonds_broken
        bonds_broken = fragment_complex.bonds_broken
        fragment_hashes = fragment_complex.fragment_hashes
        number_of_fragments = fragment_complex.number_of_fragments

        fragment_complex_ctype=molecule_entry_ctype.fragment_data[f_idx]
        fragment_complex_ctype.number_of_bonds_broken = number_of_bonds_broken
        fragment_complex_ctype.number_of_fragments = number_of_fragments

        for i, bond in enumerate(bonds_broken):
            fragment_complex_ctype.bonds_broken[i][0] = bond[0]
            fragment_complex_ctype.bonds_broken[i][1] = bond[1]

        for i, s in enumerate(fragment_hashes):
            fragment_complex_ctype.fragment_hashes[i] = s.encode('utf-8')

        molecule_entry_ctype.fragment_data[f_idx]= fragment_complex_ctype
        number_of_fragment_data+=1

    molecule_entry_ctype.number_of_fragment_data = number_of_fragment_data

    return molecule_entry_ctype




MAX_LIST_SIZE = 29 #和cpp文件同步修改
#MAX_HASH_STR_SIZE = 100

class FragmentComplex_c_type(Structure):
    _fields_ = [
                ("number_of_bonds_broken", c_int),
                ("bonds_broken", c_int*2*MAX_LIST_SIZE),
                ("number_of_fragments", c_int),
                ("fragment_hashes", c_char_p*MAX_LIST_SIZE)]

class MoleculeEntry_c_type(Structure):
    _fields_ = [
        ("number_of_fragment_data",c_int),
        ("fragment_data", FragmentComplex_c_type*MAX_LIST_SIZE),
    ]


class Return_c_type(Structure):
    _fields_ = [
        ("r", ctypes.c_bool),
        ("reactant_fragment_count", c_int),
        ("product_fragment_count", c_int),
        ("num_reactant_bonds_broken", c_int),
        ("num_product_bonds_broken", c_int),
        ("reactant_bonds_broken", c_int*2*2*MAX_LIST_SIZE),
        ("product_bonds_broken", c_int*2*2*MAX_LIST_SIZE),
        ("hashes_key", c_char_p*MAX_LIST_SIZE),
        ("hashes_value", c_int*MAX_LIST_SIZE),
        ("num_hashes", c_int)]

def cpp_function( reaction, mol_entries,lib):

    number_of_reactants=reaction['number_of_reactants']
    number_of_products=reaction['number_of_products']
    reactant0_id=reaction['reactants'][0]
    reactant1_id=reaction['reactants'][1]
    product0_id=reaction['products'][0]
    product1_id=reaction['products'][1]

    reactant0_mol_entry_ctype=create_molecule_entry(mol_entries,reactant0_id)
    reactant1_mol_entry_ctype=create_molecule_entry(mol_entries,reactant1_id)
    product0_mol_entry_ctype=create_molecule_entry(mol_entries,product0_id)
    product1_mol_entry_ctype=create_molecule_entry(mol_entries,product1_id)

    if reactant0_mol_entry_ctype is None:
        print(f"skip:mol_id:{reactant0_id}")
        return
    if reactant1_mol_entry_ctype is None:
        print(f"skip:mol_id:{reactant1_id}")
        return
    if product0_mol_entry_ctype is None:
        print(f"skip:mol_id:{product0_id}")
        return
    if product1_mol_entry_ctype is None:
        print(f"skip:mol_id:{product1_id}")
        return


    t_time=time()
    res=lib.fragment_matching_found(number_of_reactants,number_of_products,
                                ctypes.pointer(reactant0_mol_entry_ctype),
                                ctypes.pointer(reactant1_mol_entry_ctype),
                                ctypes.pointer(product0_mol_entry_ctype),
                                ctypes.pointer(product1_mol_entry_ctype)
                                )
    spend=time() - t_time
    return res.r,spend
    #print(f"res.r:{res.r},cpp_time:{time() - t_time}")

def ori_function( reaction, mols):

        reactant_fragment_indices_list = []
        product_fragment_indices_list = []

        if reaction['number_of_reactants'] == 1:
            reactant = mols[reaction['reactants'][0]]
            for i in range(len(reactant.fragment_data)):
                #print("number_of_reactants:" + str(reaction['number_of_reactants']) + " " + str(i))
                reactant_fragment_indices_list.append([i])
            # print(f"number_of_reactants:{reaction['number_of_reactants']},len(reactant_fragment_indices_list):{len(reactant_fragment_indices_list)}")


        if reaction['number_of_reactants'] == 2:
            reactant_0 = mols[reaction['reactants'][0]]
            reactant_1 = mols[reaction['reactants'][1]]
            for i in range(len(reactant_0.fragment_data)):
                for j in range(len(reactant_1.fragment_data)):
                    if (reactant_0.fragment_data[i].number_of_bonds_broken +
                        reactant_1.fragment_data[j].number_of_bonds_broken <= 1):

                        reactant_fragment_indices_list.append([i,j])
                        #print("number_of_reactants:" + str(reaction['number_of_reactants']) + " " + str(i) + " " + str(j))

            # print(f"number_of_reactants:{reaction['number_of_reactants']},len(reactant_fragment_indices_list):{len(reactant_fragment_indices_list)}")

        if reaction['number_of_products'] == 1:
            product = mols[reaction['products'][0]]
            for i in range(len(product.fragment_data)):
                #print("number_of_products:" + str(reaction['number_of_products']) + " " + str(i))
                product_fragment_indices_list.append([i])

            # print(f"number_of_products:{reaction['number_of_products']},len(product_fragment_indices_list):{len(product_fragment_indices_list)}")

        if reaction['number_of_products'] == 2:
            product_0 = mols[reaction['products'][0]]
            product_1 = mols[reaction['products'][1]]
            for i in range(len(product_0.fragment_data)):
                for j in range(len(product_1.fragment_data)):
                    if (product_0.fragment_data[i].number_of_bonds_broken +
                        product_1.fragment_data[j].number_of_bonds_broken <= 1):

                        product_fragment_indices_list.append([i,j])
                        #print("number_of_products:" + str(reaction['number_of_products']) + " " + str(i) + " " + str(j))

            # print(f"number_of_products:{reaction['number_of_products']},len(product_fragment_indices_list):{len(product_fragment_indices_list)}")

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
                        reactant_bonds_broken.append(
                            [(reactant_index, x) for x in bond])

                    for i in range(fragment_complex.number_of_fragments):
                        reactant_fragment_count += 1
                        tag = fragment_complex.fragment_hashes[i]
                        #print("reactant_hashes:" + str(reactant_fragment_indices) + " " + str(product_fragment_indices) + " " + str(frag_complex_index) + " " + tag)

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
                        product_bonds_broken.append(
                            [(product_index, x) for x in bond])


                    for i in range(fragment_complex.number_of_fragments):
                        product_fragment_count += 1
                        tag = fragment_complex.fragment_hashes[i]
                        #print("product_hashes:" + str(reactant_fragment_indices) + " " + str(product_fragment_indices) + " " + str(frag_complex_index) + " " + tag)
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
                # for i,k in enumerate(reactant_hashes):
                #     print(f"reactant_hashes:{i}" + str(k) + " " + str(reactant_hashes[k]))
                # for i,k in enumerate(product_hashes):
                #     print(f"product_hashes:{i}" + str(k) + " " + str(product_hashes[k]))
                # print(f"reactant_hashes == product_hashes:{reactant_hashes == product_hashes}")


                if reactant_hashes == product_hashes:
                    reaction['reactant_bonds_broken'] = reactant_bonds_broken
                    reaction['product_bonds_broken'] = product_bonds_broken
                    reaction['hashes'] = reactant_hashes
                    reaction['reactant_fragment_count'] = reactant_fragment_count
                    reaction['product_fragment_count'] = product_fragment_count

                    return True

        return False

def main():
    lib = ctypes.cdll.LoadLibrary("/root/HiPRGen/HiPRGen/fragment_matching_found.so")
    #定义函数参数类型和返回值类型
    lib.fragment_matching_found.argtypes = [ctypes.c_int, ctypes.c_int,
                        ctypes.POINTER(MoleculeEntry_c_type),
                        ctypes.POINTER(MoleculeEntry_c_type),
                        ctypes.POINTER(MoleculeEntry_c_type),
                        ctypes.POINTER(MoleculeEntry_c_type)]
    lib.fragment_matching_found.restype = Return_c_type


    with open("old_lib/mol_entries.pickle", 'rb') as f:
        mol_entries = pickle.load(f)

    max_list_size = 0
    for i in range(len(mol_entries)):
        mol_entry = mol_entries[i]
        max_list_size = max(max_list_size, len(mol_entry.fragment_data))
        for f_idx, fragment_complex in enumerate(mol_entry.fragment_data):
            number_of_bonds_broken = fragment_complex.number_of_bonds_broken
            number_of_fragments = fragment_complex.number_of_fragments
            max_list_size = max(max_list_size, number_of_bonds_broken, number_of_fragments)
    print("max_list_size:",max_list_size)


    for i in range(len(mol_entries)):
        for j in range(len(mol_entries)):
            reactant0_id=i
            reactant1_id=-1
            product0_id=j
            product1_id=-1

            if reactant1_id==-1:
                number_of_reactants=1
            else:
                number_of_reactants=2
            if product1_id==-1:
                number_of_products=1
            else:
                number_of_products=2
            reaction={
                'number_of_reactants':number_of_reactants,
                'number_of_products':number_of_products,
                'reactants':[reactant0_id,reactant1_id],
                'products':[product0_id,product1_id],

            }

            t_time=time()
            py_res=ori_function(copy.deepcopy(reaction),mol_entries)
            py_spend=time()-t_time

            cpp_res,cpp_spend=cpp_function(reaction,mol_entries,lib)
            if cpp_res==py_res:
                print(f"{cpp_res} py_spend:{py_spend},cpp_spend:{cpp_spend}")
            # print(f"[reactant0_id,reactant1_id]:{[reactant0_id, reactant1_id]}")
            # print(f"[product0_id,product1_id]:{[product0_id, product1_id]}")
            if cpp_res!=py_res:
                print(f"{py_res} {cpp_res} py_spend:{py_spend},cpp_spend:{cpp_spend} [reactant0_id,reactant1_id]:{[reactant0_id, reactant1_id]}")
            #print("*" * 100)


if __name__ == '__main__':
    main()