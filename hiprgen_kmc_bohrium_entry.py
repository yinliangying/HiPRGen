import os
import shutil
from pathlib import Path

from HiPRGen.launching_entry import run_with_id, run_with_smiles
from dp.launching.cli import to_runner, SubParser, default_minimal_exception_handler, run_sp_and_exit
from dp.launching.typing import BaseModel, Field
from dp.launching.typing import InputFilePath, OutputDirectory, MinMaxRange
from dp.launching.typing import Int, Float, List, Enum, String, Set, Union, Literal, Optional
from dp.launching.typing.io import InputMoleculeContent


class Subsidiary_Molecule_Options(String, Enum):
    """
    Define the target to be predicted.
    """
    Li = 'Li+1'
    OH = 'OH-1'
    C2H4 = 'C2H4'
    H2O = 'H2O'
    H2 = 'H2'
    CO = 'CO'
    CO2 = 'CO2'
    LiF = 'LiF'
    nll = 'None'
    LiCO3 = 'LiCO3-1'
    H = 'H+1'
    Li2CO3 = 'Li2CO3+1'


class SMILES(BaseModel):
    type: Literal['SMILES']
    main_mol: Optional[InputMoleculeContent] = Field(default='CCOC(=O)OCC', content_type="smiles", description="Input primary molecule")
    sub_mol: Subsidiary_Molecule_Options = Field(..., description='Input subsidiary molecule')
    n_sim: Int = Field(default=1000, description='Number of simulations for kMC')


def SMILES_task(opts: SMILES, output_directory):
    run_with_smiles(
        simulation_times=opts.n_sim,
        num_cores=8,
        main_mol_smiles=opts.main_mol.get_value(),
        sub_mol_name=opts.sub_mol.value,
        output_dir=output_directory
    )


class ID(BaseModel):
    type: Literal['ID']
    main_mol_id: Int = Field(default=10442, description='Input ID for the main molecule')
    sub_mol_id: Int = Field(default=5253, description='Input ID for the subsidiary molecule')
    n_sim: Int = Field(default=1000, description='Number of simulations for kMC')


def ID_task(opts: ID, output_directory):
    run_with_id(
        main_mol_id=opts.main_mol_id,
        sub_mol_id=opts.sub_mol_id,
        simulation_times=opts.n_sim,
        num_cores=8,
        output_dir=output_directory
    )


class kMC_pathfinding_Model(BaseModel):
    Input_Format: Union[ID, SMILES] = Field(discriminator="type")
    output_directory: OutputDirectory = Field(default='./output')


def kMC_pathfinding_task(opts: kMC_pathfinding_Model):
    if opts.Input_Format.type == 'ID':
        ID_task(opts.Input_Format, opts.output_directory.get_full_path())
    elif opts.Input_Format.type == 'SMILES':
        SMILES_task(opts.Input_Format, opts.output_directory.get_full_path())
    else:
        raise NotImplementedError


def to_parser():
    return {
        "kMC_pathfinding": SubParser(kMC_pathfinding_Model, kMC_pathfinding_task,
                                     "Run kMC and generate pathways to possiable products."),
    }


if __name__ == '__main__':
    import sys

    run_sp_and_exit(
        to_parser(),
        description="Example",
        version="0.1.0",
        exception_handler=default_minimal_exception_handler,
    )
