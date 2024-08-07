
//const int MAX_LIST_SIZE = 29; 由python传入并编译
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

// g++ -shared  -O3  -fPIC fragment_matching_found.cpp -o /root/HiPRGen/HiPRGen/fragment_matching_found.so

#include <cstring> // For strcpy()
#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <chrono>

using namespace std;
typedef struct {
    int number_of_bonds_broken;
    int (*bonds_broken)[2];
    int number_of_fragments;
    char ** fragment_hashes ;
} FragmentComplex;

typedef struct {
    int number_of_fragment_data;
    FragmentComplex *fragment_data ;
} MoleculeEntry;

typedef struct {
    bool r;
    int reactant_fragment_count;
    int product_fragment_count;
    int num_reactant_bonds_broken;
    int num_product_bonds_broken;
    int (*reactant_bonds_broken)[2][2];
    int (*product_bonds_broken)[2][2];
    char** hashes_key;
    int *hashes_value;
    int num_hashes;
} Return;

struct CharPtrHash {
    std::size_t operator()(const char* str) const {
        return std::hash<std::string>()(str);
    }
};

struct CharPtrEqual {
    bool operator()(const char* lhs, const char* rhs) const {
        return std::strcmp(lhs, rhs) == 0;
    }
};

bool areMapsEqual(const std::unordered_map< char*, int, CharPtrHash, CharPtrEqual> & map1, const std::unordered_map< char*, int, CharPtrHash, CharPtrEqual> & map2) {
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
    auto start = std::chrono::high_resolution_clock::now();


    vector<pair<int, int>> reactant_fragment_indices_list;
    vector<pair<int, int>> product_fragment_indices_list;
    Return result;
    result.r = false;

    if (number_of_reactants == 1) {
        for (int i = 0; i < reactant0_mol->number_of_fragment_data; i++) {
        }
    }
    else if (number_of_reactants == 2) {
        for (int i = 0; i < reactant0_mol->number_of_fragment_data; i++) {
            for (int j = 0; j < reactant1_mol->number_of_fragment_data; j++) {
                if ((reactant0_mol->fragment_data[i].number_of_bonds_broken +
                        reactant1_mol->fragment_data[j].number_of_bonds_broken) <= 1) {
                    reactant_fragment_indices_list.push_back(make_pair(i, j));
                }
            }
        }
    }

    if (number_of_products == 1) {
        for (int i = 0; i < product0_mol->number_of_fragment_data; i++) {
            product_fragment_indices_list.push_back(make_pair(i, -1));
        }
    }
    else if (number_of_products == 2) {
        for (int i = 0; i < product0_mol->number_of_fragment_data; i++) {
            for (int j = 0; j < product1_mol->number_of_fragment_data; j++) {
                if ((product0_mol->fragment_data[i].number_of_bonds_broken +
                        product1_mol->fragment_data[j].number_of_bonds_broken) <= 1) {
                    product_fragment_indices_list.push_back(make_pair(i, j));
                }
            }
        }
    }

    for (int tmp_reactant_fragment_idx = 0; tmp_reactant_fragment_idx < reactant_fragment_indices_list.size(); tmp_reactant_fragment_idx++){
        for (int tmp_product_fragment_idx = 0; tmp_product_fragment_idx < product_fragment_indices_list.size(); tmp_product_fragment_idx++) {
            int reactant_fragment_count = 0;
            int product_fragment_count = 0;
//            int  reactant_bonds_broken[MAX_LIST_SIZE][2][2];
//            int reactant_bonds_broken_len=0;
//            int  product_bonds_broken[MAX_LIST_SIZE][2][2];
//            int product_bonds_broken_len=0;
            std::unordered_map< char*, int, CharPtrHash, CharPtrEqual>  reactant_hashes;
            std::unordered_map< char*, int, CharPtrHash, CharPtrEqual>  product_hashes;

            int product_fragment_indices[2]={product_fragment_indices_list[tmp_product_fragment_idx].first,product_fragment_indices_list[tmp_product_fragment_idx].second};
            for (int product_index = 0; product_index < 2; product_index++) {

                int frag_complex_index = product_fragment_indices[product_index];

                FragmentComplex *fragment_complex;
                if (product_index == 0) {
                    fragment_complex = &(product0_mol->fragment_data[frag_complex_index]);
                }
                else if (product_index == 1) {
                    if (frag_complex_index == -1){
                        continue;
                    }
                    fragment_complex = &(product1_mol->fragment_data[frag_complex_index]);
                }
                else {
                    continue;
                }

//                for (int i = 0; i < fragment_complex->number_of_bonds_broken;i++) {
//                    product_bonds_broken[product_bonds_broken_len][0][0]=product_index;
//                    product_bonds_broken[product_bonds_broken_len][0][1]=fragment_complex->bonds_broken[i][0];
//                    product_bonds_broken[product_bonds_broken_len][1][0]=product_index;
//                    product_bonds_broken[product_bonds_broken_len][1][1]=fragment_complex->bonds_broken[i][1];
//                    product_bonds_broken_len++;
//                }

                for (int i = 0; i < fragment_complex->number_of_fragments; i++) {
                    product_fragment_count++;
                    char* tag = fragment_complex->fragment_hashes[i];
                    if (product_hashes.count(tag) > 0) {
                        product_hashes[tag]++;
                    }
                    else {
                        product_hashes[tag] = 1;
                    }
                }
            }


            int reactant_fragment_indices[2]={reactant_fragment_indices_list[tmp_reactant_fragment_idx].first,reactant_fragment_indices_list[tmp_reactant_fragment_idx].second};
            for (int reactant_index = 0; reactant_index < 2; reactant_index++) {
                int frag_complex_index = reactant_fragment_indices[reactant_index];

                FragmentComplex *fragment_complex;
                if (reactant_index == 0) {
                    fragment_complex =&( reactant0_mol->fragment_data[frag_complex_index]);
                }
                else if (reactant_index == 1) {
                    if (frag_complex_index == -1){
                        continue;
                    }
                    fragment_complex = &(reactant1_mol->fragment_data[frag_complex_index]);
                }
                else {
                    continue;
                }

//                for (int i = 0; i < fragment_complex->number_of_bonds_broken;i++)  {
//                    reactant_bonds_broken[reactant_bonds_broken_len][0][0]=reactant_index;
//                    reactant_bonds_broken[reactant_bonds_broken_len][0][1]=fragment_complex->bonds_broken[i][0];
//                    reactant_bonds_broken[reactant_bonds_broken_len][1][0]=reactant_index;
//                    reactant_bonds_broken[reactant_bonds_broken_len][1][1]=fragment_complex->bonds_broken[i][1];
//                    reactant_bonds_broken_len++;
//                }

                for (int i = 0; i < fragment_complex->number_of_fragments; i++) {
                    reactant_fragment_count++;
                    char* tag = fragment_complex->fragment_hashes[i];

                    if (reactant_hashes.count(tag) > 0) {
                        reactant_hashes[tag]++;
                    }
                    else {
                        reactant_hashes[tag] = 1;
                    }
                }
            }



            if (number_of_reactants == 2 && number_of_products == 2 &&
                reactant_fragment_count == 2 && product_fragment_count == 2) {
                continue;
            }

            bool isEqual = areMapsEqual(reactant_hashes, product_hashes);

            if (isEqual) {
                result.r = true;
                return result;  // 只返回是否通过，不返回其他信息了，因为这个函数通过率非常低，如果返回true让python计算其他结果就可以
//                result.reactant_fragment_count = reactant_fragment_count;
//                result.product_fragment_count = product_fragment_count;
//                result.num_reactant_bonds_broken = reactant_bonds_broken_len;
//                result.num_product_bonds_broken =product_bonds_broken_len;
//                result.num_hashes = reactant_hashes.size();
//
//
//                for (int i=0; i<reactant_bonds_broken_len; i++){
//
//                    result.reactant_bonds_broken[i][0][0] =reactant_bonds_broken[i][0][0];
//                    result.reactant_bonds_broken[i][0][1] = reactant_bonds_broken[i][0][1];
//                    result.reactant_bonds_broken[i][1][0] = reactant_bonds_broken[i][1][0];
//                    result.reactant_bonds_broken[i][1][1] = reactant_bonds_broken[i][1][1];
//                }
//
//                for (int i=0; i<product_bonds_broken_len; i++){
//
//
//                    result.product_bonds_broken[i][0][0] = product_bonds_broken[i][0][0];
//                    result.product_bonds_broken[i][0][1] = product_bonds_broken[i][0][1];
//                    result.product_bonds_broken[i][1][0] = product_bonds_broken[i][1][0];
//                    result.product_bonds_broken[i][1][1] = product_bonds_broken[i][1][1];
//                }
//
//                int index = 0;
//                for (auto& hash : reactant_hashes) {
//                    result.hashes_key[index] =hash.first;
//                    result.hashes_value[index] = hash.second;
//                    index++;
//                }
//
//                return result;
            }
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    duration = duration;
    return result;
}