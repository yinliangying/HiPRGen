# import pandas as pd
# import pickle
# template_list=[]
# with open("local_mapper_templates_202403.pkl","rb") as fp:
#     template_set=pickle.load(fp)
#
# for template in template_set:
#     template_list.append(template)
#
# #https://github.com/hesther/templatecorr
# data = pd.read_hdf("uspto_460k_unique_templates.hdf5", "table")
#
# for i in range(len(data)):
#     products, reactants =data["retro_template"][i].split(">>")
#     template_str=">>".join([reactants, products])
#     template_list.append(template_str)
#
#
# with open("templatecorr_templates_202411.pkl","wb") as fp:
#     pickle.dump(template_list,fp)



from rdkit import Chem
from rdkit.Chem import AllChem

# 定义反应模板
template_srarts = '[C:1]-[O:2]>>[C:1]-[S:2]'
template_obj= AllChem.ReactionFromSmarts(template_srarts)

def is_reaction_matched(template_obj, reactant_mols, product_smiles_set):
    template_products_mols = template_obj.RunReactants(reactant_mols)[0]
    template_products_set = set([Chem.MolToSmiles(mol) for mol in template_products_mols])
    return template_products_set == product_smiles_set

# 示例数据
tmp_row = {
    'substrates': ['CCO'],
    'products': ['CCS']
}

# 判断反应是否匹配模板
reactant_mols=[Chem.MolFromSmiles(s) for s in tmp_row['substrates']]
product_smiles_set=set(tmp_row['products'])
matched = is_reaction_matched(template_obj, reactant_mols, product_smiles_set)
print(f"Reaction matched: {matched}")
