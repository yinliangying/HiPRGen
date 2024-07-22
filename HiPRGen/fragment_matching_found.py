
import ctypes
import pickle
from ctypes import Structure, c_int, POINTER,c_char,c_char_p,c_bool
from HiPRGen.mol_entry import MoleculeEntry
from time import time
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


def main():
    lib = ctypes.cdll.LoadLibrary("/root/HiPRGen/HiPRGen/fragment_matching_found.so")

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


    reactant0_id=1
    reactant1_id=2
    product0_id=3
    product1_id=-1
    number_of_reactants=2
    number_of_products=1

    reactant0_mol_entry_ctype=create_molecule_entry(mol_entries,reactant0_id)
    reactant1_mol_entry_ctype=create_molecule_entry(mol_entries,reactant1_id)
    product0_mol_entry_ctype=create_molecule_entry(mol_entries,product0_id)
    product1_mol_entry_ctype=create_molecule_entry(mol_entries,product1_id)
    if reactant0_mol_entry_ctype is None or \
            reactant1_mol_entry_ctype is None or \
            product0_mol_entry_ctype is None or \
            product1_mol_entry_ctype is None:
        print("skip")
        exit()



    #定义函数参数类型和返回值类型
    lib.fragment_matching_found.argtypes = [ctypes.c_int, ctypes.c_int,
                        ctypes.POINTER(MoleculeEntry_c_type),
                        ctypes.POINTER(MoleculeEntry_c_type),
                        ctypes.POINTER(MoleculeEntry_c_type),
                        ctypes.POINTER(MoleculeEntry_c_type)]

    lib.fragment_matching_found.restype = Return_c_type
    t_time=time()
    r=lib.fragment_matching_found(number_of_reactants,number_of_products,
                                ctypes.pointer(reactant0_mol_entry_ctype),
                                ctypes.pointer(reactant1_mol_entry_ctype),
                                ctypes.pointer(product0_mol_entry_ctype),
                                ctypes.pointer(product1_mol_entry_ctype)
                                )
    print(r.r,time()-t_time)



if __name__ == '__main__':
    main()