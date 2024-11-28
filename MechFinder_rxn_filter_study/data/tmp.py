import pandas as pd
# import pickle
# template_list=[]
# with open("local_mapper_templates_202403.pkl","rb") as fp:
#     template_set=pickle.load(fp)
#
# for template in template_set:
#     template_list.append(template)
#
#https://github.com/hesther/templatecorr
data = pd.read_hdf("uspto_460k_unique_templates.hdf5", "table")
print(len(data))
exit()
# for i in range(len(data)):
#     products, reactants =data["retro_template"][i].split(">>")
#     template_str=">>".join([reactants, products])
#     template_list.append(template_str)
#
#
# with open("templatecorr_templates_202411.pkl","wb") as fp:
#     pickle.dump(template_list,fp)

from rdkit import Chem
from rdkit.Chem import Draw
mol = Chem.MolFromSmiles('[CH2][C]OC[CH-]OP(O)(F)(F)F')

# 遍历分子中的原子并移除自由基
for atom in mol.GetAtoms():
    if atom.GetNumRadicalElectrons() > 0:
        atom.SetNumRadicalElectrons(0)

# 规范化 SMILES 字符串
smiles = Chem.MolToSmiles(mol, canonical=True)
print(smiles)
mol=Chem.MolFromSmiles("[CH2]COC[CH-]OP(O)(F)(F)F")
mol=Chem.RemoveHs(mol)
for atom in mol.GetAtoms():
    if atom.GetNumRadicalElectrons() > 0:
        atom.SetNumRadicalElectrons(0)
    atom.SetNumRadicalElectrons(0)
mol=Chem.AddHs(mol)
img = Draw.MolToImage(mol)
print(Chem.MolToSmiles(mol))
mol=Chem.RemoveHs(mol)
print(Chem.MolToSmiles(mol))
img.show()


from rdkit import Chem
from rdkit.Chem import AllChem
import pickle

with open("templatecorr_templates_202411.pkl","rb") as fp:
    template_smarts_list=pickle.load(fp)
template_smarts_obj_list = []
for template_str in template_smarts_list:
    template_obj = AllChem.ReactionFromSmarts(template_str)
    template_smarts_obj_list.append((template_str, template_obj))

def is_reaction_matched(template_obj, reactant_mols, product_smiles_set):
    template_products_mols = template_obj.RunReactants(reactant_mols)[0]
    template_products_set = set([Chem.MolToSmiles(mol) for mol in template_products_mols])
    print(".".join([AllChem.MolToSmiles(mol) for mol in reactant_mols]+list(template_products_set)+list(product_smiles_set)))

    return template_products_set == product_smiles_set



input_str="""[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , OP(F)(F)(F)O[CH]CO[C]=[C]F.C=COC[CH-]OP(O)(F)(F)F>>[CH2][C]OC[CH-]OP(O)(F)(F)F.OP(F)(F)(F)O[CH]CO[CH][C]F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , OP(F)(F)(F)OC=CO[CH-]CF.[C]=COC[CH]OP(O)(F)(F)F>>OP(F)(F)(F)O[CH-][CH]O[C]CF.[CH][CH]OC[CH]OP(O)(F)(F)F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , OP(F)(F)(F)OC=CO[CH-]CF.[C]=COC[CH]OP(O)(F)(F)F>>OP(F)(F)(F)O[CH-][C]O[C]CF.[CH][CH]OC[CH]OP(O)(F)(F)F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , C=COC[CH-]OP(O)(F)(F)F.OP(F)(F)(F)OC=[C]O[C]CF>>OP(F)(F)(F)O[CH][C]O[CH]CF.[CH2][C]OC[CH-]OP(O)(F)(F)F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , O=P(F)(F)OC=[C]OCCF.OP(F)(F)(F)O[C-]=COC[CH]F>>O=P(F)(F)OC[C]OCCF.OP(F)(F)(F)O[CH-][CH]OC[C]F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , [CH]=COC=COP(O)(F)(F)F.OP(F)(F)(F)O[CH-]CO[C]CF>>[O-]P(F)(F)(F)O[C]CO[C]CF.[CH][CH]OC[CH]OP(O)(F)(F)F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , [CH]=COC=COP(O)(F)(F)F.OP(F)(F)(F)O[CH-]COC[C]F>>[O-]P(F)(F)(F)O[C]COC[C]F.[CH][CH]OC[CH]OP(O)(F)(F)F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , C=COC=COP(O)(F)(F)F.OP(F)(F)(F)O[CH]CO[C]=[C]F>>[CH2][CH+]OC=COP(O)(F)(F)F.OP(F)(F)(F)O[CH]CO[CH][C]F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , C=COC=COP(O)(F)(F)F.OP(F)(F)(F)OC=[C]O[C]CF>>[CH2][CH+]OC=COP(O)(F)(F)F.OP(F)(F)(F)O[CH][C]O[CH]CF
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , C=COC=COP(O)(F)(F)F.OP(F)(F)(F)O[C]=CO[C]CF>>[CH2][CH+]OC=COP(O)(F)(F)F.OP(F)(F)(F)O[C][CH]O[CH]CF
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , F[C]CO[CH-]OP(F)F.OP(F)(F)(F)OC=[C]OC=CF>>F[C][CH]O[CH-]OP(F)F.OP(F)(F)(F)O[CH][C]O[CH]CF
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , F[C]CO[CH-]OP(F)F.OP(F)(F)(F)O[C]=COC=CF>>OP(F)(F)(F)O[C][CH]O[CH]CF.F[C][CH]O[CH-]OP(F)F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , C=COC=[C]OP(O)(F)(F)F.[O]P(F)(F)(F)OCCO[C]CF>>OP(F)(F)(F)O[C][CH]O[CH]CF.[CH2][C]OCCOP([O])(F)(F)F
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , FC[C]O[CH-]OP(F)F.OP(F)(F)(F)OC=[C]OC=CF>>F[CH][C]O[CH-]OP(F)F.OP(F)(F)(F)O[CH][C]O[CH]CF
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , FC[C]O[CH-]OP(F)F.OP(F)(F)(F)O[C]=COC=CF>>F[CH][C]O[CH-]OP(F)F.OP(F)(F)(F)O[C][CH]O[CH]CF
[C:1]=[C:2].[C:3]=[C:4]>>[C:1]-[C:2].[C:3]-[C:4] , F[C]=CO[C]OP(F)F.OP(F)(F)(F)O[C-]=COC[CH]F>>F[C]CO[C]OP(F)F.OP(F)(F)(F)O[CH-][CH]OC[C]F"""
for line in input_str.split("\n"):
    result_template_smarts, rxn_smiles = line.split(" , ")

    reactant_smiles_list=[]
    product_smiles_list=[]
    for smiles in rxn_smiles.split(">>")[0].split("."):
        mol = AllChem.MolFromSmiles(smiles)
        for idx, atom in enumerate(mol.GetAtoms()):
            mol.GetAtomWithIdx(idx).SetNumRadicalElectrons(0)
        no_radical_smiles = Chem.MolToSmiles(mol)
        reactant_smiles_list.append(no_radical_smiles)
    for smiles in rxn_smiles.split(">>")[1].split("."):
        mol = AllChem.MolFromSmiles(smiles)
        for idx, atom in enumerate(mol.GetAtoms()):

            atom.SetNumRadicalElectrons(0)
        no_radical_smiles = Chem.MolToSmiles(mol)
        product_smiles_list.append(no_radical_smiles)

    reactant_mols = [Chem.MolFromSmiles(smiles) for smiles in reactant_smiles_list] # [Chem.MolFromSmiles(smiles) for smiles in rxn_smiles.split(">>")[0].split(".")]
    product_smiles_set =  set(product_smiles_list)   # set(rxn_smiles.split(">>")[1].split("."))

    matched = False
    matched_template_smarts = ""
    for template_smarts, template_obj in template_smarts_obj_list:
        if template_smarts != result_template_smarts:
            continue

        template_reactant_num = len(template_smarts.split(">>")[0].split("."))
        template_product_num = len(template_smarts.split(">>")[1].split("."))
        if len(rxn_smiles.split(">>")[0].split(".")) == template_reactant_num and len(
                rxn_smiles.split(">>")[1].split(".")) == template_product_num:
            if is_reaction_matched(template_obj, reactant_mols, product_smiles_set):
                matched = True
                matched_template_smarts = template_smarts
                break
            print()
    if matched:
        print()
