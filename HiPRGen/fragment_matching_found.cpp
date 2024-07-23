// python 调用方式
//number_of_reactants=1,
//number_of_products=2
//
//reactant0_mol=class MoleculeEntry{fragment_data:[class FragmentComplex]}
//
//class MoleculeEntry
//{
//    list fragment_data =[
//    class FragmentComplex{
//        int number_of_bonds_broken
//        list bonds_broken=[(0,1),...]
//        int number_of_fragments
//        fragment_hashes=[str,str,...]
//        }
//    ]
//
//}

// g++ -shared -fPIC fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so

#include <cstring> // For strcpy()
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>

const int MAX_LIST_SIZE = 29;

typedef struct {
    int number_of_bonds_broken;
    int bonds_broken[MAX_LIST_SIZE][2];
    int number_of_fragments;
    char* fragment_hashes[MAX_LIST_SIZE];
} FragmentComplex;

typedef struct {
    int number_of_fragment_data;
    FragmentComplex fragment_data[MAX_LIST_SIZE];
} MoleculeEntry;

typedef struct {
    bool r;
    int reactant_fragment_count;
    int product_fragment_count;
    int num_reactant_bonds_broken;
    int num_product_bonds_broken;
    int reactant_bonds_broken[MAX_LIST_SIZE][2][2];
    int product_bonds_broken[MAX_LIST_SIZE][2][2];
    char* hashes_key[MAX_LIST_SIZE];
    int hashes_value[MAX_LIST_SIZE];
    int num_hashes;
} Return;


bool areMapsEqual(const std::unordered_map<std::string, int>& map1, const std::unordered_map<std::string, int>& map2) {
    if (map1.size() != map2.size()) {
        return false;
    }

    for (const auto& pair : map1) {
        auto it = map2.find(pair.first);
        if (it == map2.end() || it->second != pair.second) {
            return false;
        }
    }

    return true;
}

extern "C" Return fragment_matching_found(int number_of_reactants, int number_of_products,
                MoleculeEntry *reactant0_mol, MoleculeEntry *reactant1_mol,
                MoleculeEntry *product0_mol, MoleculeEntry *product1_mol) {
    std::vector<std::vector<int>> reactant_fragment_indices_list;
    std::vector<std::vector<int>> product_fragment_indices_list;
    Return result;

    if (number_of_reactants == 1) {
        for (int i = 0; i < reactant0_mol->number_of_fragment_data; i++) {
//            std::cout << "number_of_reactants:" << number_of_reactants<<" "<<i << std::endl;
            reactant_fragment_indices_list.push_back({i, -1});
        }
        std::cout << "number_of_reactants:" << number_of_reactants<<" "<<reactant_fragment_indices_list.size() << std::endl;
    }

    if (number_of_reactants == 2) {
        for (int i = 0; i < reactant0_mol->number_of_fragment_data; i++) {
            for (int j = 0; j < reactant1_mol->number_of_fragment_data; j++) {
                if ((reactant0_mol->fragment_data[i].number_of_bonds_broken +
                        reactant1_mol->fragment_data[j].number_of_bonds_broken) <= 1) {
                    reactant_fragment_indices_list.push_back({i, j});
//                    std::cout << "number_of_reactants:" << number_of_reactants<<" "<<i<<j << std::endl;
                }
            }
        }
        std::cout << "number_of_reactants:" << number_of_reactants<<" "<<reactant_fragment_indices_list.size() << std::endl;
    }

    if (number_of_products == 1) {
        for (int i = 0; i < product0_mol->number_of_fragment_data; i++) {
//            std::cout << "number_of_products:" << number_of_products<<" "<<i << std::endl;
            product_fragment_indices_list.push_back({i, -1});
        }
        std::cout << "number_of_products:" << number_of_products<<" "<<product_fragment_indices_list.size() << std::endl;

    }

    if (number_of_products == 2) {
        for (int i = 0; i < product0_mol->number_of_fragment_data; i++) {
            for (int j = 0; j < product1_mol->number_of_fragment_data; j++) {
                if ((product0_mol->fragment_data[i].number_of_bonds_broken +
                        product1_mol->fragment_data[j].number_of_bonds_broken) <= 1) {
                    product_fragment_indices_list.push_back({i, j});
//                    std::cout << "number_of_products:" << number_of_products<<" "<<i<<j << std::endl;
                }
            }
        }
        std::cout << "number_of_products:" << number_of_products<<" "<<product_fragment_indices_list.size() << std::endl;
    }

    for (auto& reactant_fragment_indices : reactant_fragment_indices_list) {
        for (auto& product_fragment_indices : product_fragment_indices_list) {
            int reactant_fragment_count = 0;
            int product_fragment_count = 0;
            std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> reactant_bonds_broken;
            std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> product_bonds_broken;
            std::unordered_map<std::string, int> reactant_hashes;
            std::unordered_map<std::string, int> product_hashes;

            for (int reactant_index = 0; reactant_index < reactant_fragment_indices.size(); reactant_index++) {
                int frag_complex_index = reactant_fragment_indices[reactant_index];

                FragmentComplex fragment_complex;
                if (reactant_index == 0) {
                    fragment_complex = reactant0_mol->fragment_data[frag_complex_index];
                }
                else if (reactant_index == 1) {
                    if (frag_complex_index == -1){
                        continue;
                    }
                    fragment_complex = reactant1_mol->fragment_data[frag_complex_index];
                }
                else {
                    continue;
                }

                for (auto& bond : fragment_complex.bonds_broken) {
                    reactant_bonds_broken.push_back(std::make_pair(std::make_pair(reactant_index, bond[0]), std::make_pair(reactant_index, bond[1])));
                }

                for (int i = 0; i < fragment_complex.number_of_fragments; i++) {
                    reactant_fragment_count++;
                    std::string tag = fragment_complex.fragment_hashes[i];

//                    std::cout << "reactant_hashes: ";
//                    std::cout << "[" ;
//                    for (const auto& index : reactant_fragment_indices) {
//                        std::cout << index << " ";
//                    }
//                    std::cout << "][";
//                    for (const auto& index : product_fragment_indices) {
//                        std::cout << index << " ";
//                    }
//                    std::cout <<"] " << frag_complex_index << " "  << tag << std::endl;

                    if (reactant_hashes.count(tag) > 0) {
                        reactant_hashes[tag]++;
                    }
                    else {
                        reactant_hashes[tag] = 1;
                    }
                }
            }

            for (int product_index = 0; product_index < product_fragment_indices.size(); product_index++) {
                int frag_complex_index = product_fragment_indices[product_index];

                FragmentComplex fragment_complex;
                if (product_index == 0) {
                    fragment_complex = product0_mol->fragment_data[frag_complex_index];
                }
                else if (product_index == 1) {
                    if (frag_complex_index == -1){
                        continue;
                    }
                    fragment_complex = product1_mol->fragment_data[frag_complex_index];
                }
                else {
                    continue;
                }

                for (auto& bond : fragment_complex.bonds_broken) {
                    product_bonds_broken.push_back(std::make_pair(std::make_pair(product_index, bond[0]), std::make_pair(product_index, bond[1])));
                }

                for (int i = 0; i < fragment_complex.number_of_fragments; i++) {
                    product_fragment_count++;
                    std::string tag = fragment_complex.fragment_hashes[i];

//                    std::cout << "product_hashes: ";
//                    std::cout << "[" ;
//                    for (const auto& index : reactant_fragment_indices) {
//                        std::cout << index << " ";
//                    }
//                    std::cout << "][";
//                    for (const auto& index : product_fragment_indices) {
//                        std::cout << index << " ";
//                    }
//                    std::cout <<"] " << frag_complex_index << " "  << tag << std::endl;
                    if (product_hashes.count(tag) > 0) {
                        product_hashes[tag]++;
                    }
                    else {
                        product_hashes[tag] = 1;
                    }
                }
            }


            if (number_of_reactants == 2 && number_of_products == 2 &&
                reactant_fragment_count == 2 && product_fragment_count == 2) {
                continue;
            }

            bool isEqual = areMapsEqual(reactant_hashes, product_hashes);
            if (isEqual) {
                int tmp_reactant_hashes_index = 0;
                int tmp_product_hashes_index = 0;
                for (auto& hash : reactant_hashes) {
                        std::cout<<tmp_reactant_hashes_index<< hash.first << hash.second << std::endl;
                        tmp_reactant_hashes_index++;
                    }
                std::cout<<"reactant_hashes.size()"<<reactant_hashes.size() << std::endl;
                for (auto& hash : product_hashes) {
                        std::cout<<product_hashes<< hash.first << hash.second << std::endl;
                        tmp_product_hashes_index++;
                    }
                std::cout<<"product_hashes.size() "<<product_hashes.size() << std::endl;

                std::cout<<"isEqual"<<isEqual << std::endl;
            }

            if (isEqual) {
                result.r = true;
                result.reactant_fragment_count = reactant_fragment_count;
                result.product_fragment_count = product_fragment_count;
                result.num_reactant_bonds_broken = reactant_bonds_broken.size();
                result.num_product_bonds_broken = product_bonds_broken.size();
                result.num_hashes = reactant_hashes.size();

                int index = 0;
                for (auto& bond : reactant_bonds_broken) {
                    result.reactant_bonds_broken[index][0][0] = bond.first.first;
                    result.reactant_bonds_broken[index][0][1] = bond.first.second;
                    result.reactant_bonds_broken[index][1][0] = bond.second.first;
                    result.reactant_bonds_broken[index][1][1] = bond.second.second;
                    index++;
                }

                index = 0;
                for (auto& bond : product_bonds_broken) {
                    result.product_bonds_broken[index][0][0] = bond.first.first;
                    result.product_bonds_broken[index][0][1] = bond.first.second;
                    result.product_bonds_broken[index][1][0] = bond.second.first;
                    result.product_bonds_broken[index][1][1] = bond.second.second;
                    index++;
                }

                index = 0;
                for (auto& hash : reactant_hashes) {
                    result.hashes_key[index] =new char[hash.first.size() + 1];
                    strcpy(result.hashes_key[index], hash.first.c_str()); //delete[] result.hashes_key[index]; // 不要忘记释放内存
                    result.hashes_value[index] = hash.second;
                    index++;
                }

                return result;
            }
        }
    }

    return result;
}