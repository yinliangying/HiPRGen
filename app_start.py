import os
from dp.launching.cli import SubParser, default_minimal_exception_handler, run_sp_and_exit, to_runner
from dp.launching.typing import List, BaseModel, Field, OutputDirectory, InputFilePath, Int, Union, String, Literal, \
    Field, Enum
from dp.launching.typing.io import InputMoleculeContent



def SCORING_func(output_dir,
                 method_type,
                 ):
    HiPRGen_dir="/root/HiPRGen"
    network_json_path = "%s/data/libe.json"%(HiPRGen_dir)
    os.system("mkdir -p %s" % (output_dir))
    number_of_threads = os.popen("nproc").read().strip()  # os.system("nproc")
    if method_type.contact.type == "LIBE_ID_input":
        observed_molecule_libe_id_list = method_type.contact.observed_molecule_libe_id_list
        init_molecule_libe_id_list = method_type.contact.init_molecule_libe_id_list

        if observed_molecule_libe_id_list == "":
            os.system(f" . {HiPRGen_dir}/nix_env.sh && python {HiPRGen_dir}/bohrium_start.py  \
            --network_json_path  {network_json_path}  --molecule_type libe_id  \
            --number_of_threads {number_of_threads} --work_dir  {output_dir} \
             --init_molecule_libe_id_list {init_molecule_libe_id_list}  \
            ")
        else:
            os.system(f" . {HiPRGen_dir}/nix_env.sh && python {HiPRGen_dir}/bohrium_start.py  \
            --network_json_path  {network_json_path}  --molecule_type libe_id  \
            --number_of_threads {number_of_threads} --work_dir  {output_dir} \
             --init_molecule_libe_id_list {init_molecule_libe_id_list} --observed_molecule_libe_id_list {observed_molecule_libe_id_list} \
            ")

    elif method_type.contact.type == "SMILES_input":
        observed_molecule_smiles_list = method_type.contact.observed_molecule_smiles_list
        init_molecule_smiles_list = method_type.contact.init_molecule_smiles_list
        with open("init_molecule_smiles_list.txt","w") as init_molecule_smiles_list_fp:
            init_molecule_smiles_list_fp.write(init_molecule_smiles_list)
        with open("observed_molecule_smiles_list.txt","w") as observed_molecule_smiles_list_fp:
            observed_molecule_smiles_list_fp.write(observed_molecule_smiles_list)

        #old
        if observed_molecule_smiles_list == "":
            os.system(f" . {HiPRGen_dir}/nix_env.sh && python {HiPRGen_dir}/bohrium_start.py  \
            --network_json_path  {network_json_path}  --molecule_type SMILES  \
            --number_of_threads {number_of_threads} --work_dir  {output_dir} \
             --init_molecule_smiles_list_file  init_molecule_smiles_list.txt  \
            ")
        else:
            os.system(f" . {HiPRGen_dir}/nix_env.sh && python {HiPRGen_dir}/bohrium_start.py  \
            --network_json_path  {network_json_path}  --molecule_type SMILES  \
            --number_of_threads {number_of_threads} --work_dir  {output_dir} \
              --init_molecule_smiles_list_file  init_molecule_smiles_list.txt  \
              --observed_molecule_smiles_list_file observed_molecule_smiles_list.txt  ")



    elif method_type.contact.type == "find_LIBE_ID_by_formula_alphabetical":
        formula_alphabetical_list = method_type.contact.formula_alphabetical_list.replace(" ", "_")
        os.system(f" . {HiPRGen_dir}/nix_env.sh && python {HiPRGen_dir}/bohrium_start.py  \
                    --network_json_path  {network_json_path}  --molecule_type find_formula  \
                    --number_of_threads {number_of_threads} --work_dir  {output_dir} \
                    --find_libe_id_by_formula_alphabetical_list {formula_alphabetical_list} \
                    ")


class SMILES_input(BaseModel):
    type: Literal["SMILES_input"]
    init_molecule_smiles_list: InputMoleculeContent = Field(content_type="smiles",
                                                            description="""molecule(s) with SMILES format to be analysed exp O=C1OCCO1.C=C=C=C=O . If have 2 or more molecules, they should be splited by '.'""")
    observed_molecule_smiles_list: InputMoleculeContent = Field(default="", content_type="smiles",
                                                                description="molecule(s) with SMILES format to be observed,which could be empty. If have 2 or more molecules, they should be splited by '.'")


class LIBE_ID_input(BaseModel):
    type: Literal["LIBE_ID_input"]
    init_molecule_libe_id_list: String = Field(
        description="""molecule(s) with LIBE ID to be analysed exp EC:libe-115834, which are splited by ',' if have 2 or more molecules""")
    observed_molecule_libe_id_list: String = Field(
        description="molecule(s) with LIBE ID to be observed,which could be empty and are splited by ',' if have 2 or more molecules")


class find_LIBE_ID_by_formula_alphabetical(BaseModel):
    type: Literal["find_LIBE_ID_by_formula_alphabetical"]
    formula_alphabetical_list: String = Field(
        description="a list of alphabetical formula splited by ',', which is used for find the LIBE ID of molecule in dataset,exp:C1 O2,C3 H4 Li1 O3")


class UnionType(BaseModel):
    contact: Union[SMILES_input, LIBE_ID_input, find_LIBE_ID_by_formula_alphabetical] = Field(discriminator="type")


class Options(BaseModel):
    method_type: UnionType
    output_dir: OutputDirectory = Field(
        default="./outputs"
    )  # default will be override after online


class global_opt(Options, BaseModel):
    ...


def runner(opts: global_opt) -> int:
    pass
    SCORING_func(
        output_dir=opts.output_dir, method_type=opts.method_type)


def to_parser():
    return to_runner(
        global_opt,
        runner,
        version="0.1.0",
        exception_handler=default_minimal_exception_handler,
    )


if __name__ == '__main__':
    import sys

    to_parser()(sys.argv[1:])
