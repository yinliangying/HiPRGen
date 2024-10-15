import os
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from ase.db import connect
from ase import Atoms
from ase.symbols import symbols2numbers
from openbabel import openbabel


def plot_molecule_to_pdf(mol, output_filename="molecule.pdf", charge: int=0):
    """
    Plot an RDKit molecule object to a PDF file.
    """
    # Generate 2D coordinates for the molecule
    AllChem.Compute2DCoords(mol)

    # Create a matplotlib figure
    fig, ax = plt.subplots(figsize=(8, 8))

    # Draw the molecule
    img = Draw.MolToImage(mol)

    # Display the molecule in the matplotlib figure
    ax.imshow(img)
    ax.axis('off')  # Hide axes

    # Get the dimensions of the image
    img_height, img_width = img.size

    # Add charge in the corner with lines
    if charge != 0:
        charge_str = f"{'+' if charge > 0 else '-'}{abs(charge)}"
        ax.text(img_width, -0.05*img_height, charge_str, fontsize=20,
                ha='left', va='top', color='black')

        # Add corner lines
        ax.plot([img_width * 0.7, img_width], [0, 0], color='black', linewidth=2)
        ax.plot([img_width, img_width], [0, img_height * 0.3], color='black', linewidth=2)

    # Save as PDF
    plt.savefig(output_filename, format="pdf", bbox_inches="tight", dpi=600)
    plt.close(fig)  # Close the figure to free up memory

    print(f"Molecule saved as '{output_filename}'")


def plot_atom_to_pdf(atom_symbol, output_filename="molecule.pdf", charge: int=0):
    # Create a matplotlib figure
    fig, ax = plt.subplots(figsize=(8, 8))

    # Remove axes
    ax.axis('off')

    # Determine the symbol based on the charge
    if charge == 0:
        symbol = f'{atom_symbol}'
    elif charge == 1:
        symbol = atom_symbol + '$^{+}$'
    elif charge == -1:
        symbol = atom_symbol + '$^{-}$'
    else:
        raise ValueError("Charge must be -1, 0, or 1")

    # Add the text "H" in the center
    ax.text(0.5, 0.5, symbol, fontsize=50, ha='center', va='center')

    # Save as PDF
    plt.savefig(output_filename, format="pdf", bbox_inches="tight", dpi=600)
    plt.close(fig)  # Close the figure to free up memory

    print(f"Molecule saved as '{output_filename}'")


def xyz_2_db_mol(idx):
    # Initialize Open Babel conversion
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "smi")
    # Read the XYZ file into an OBMol object
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, f'{idx}.xyz')
    mol.ConnectTheDots()
    mol.PerceiveBondOrders()
    return mol


def obmol_to_rdkit_mol(obmol):
    # Create an empty RDKit molecule
    rdkit_mol = Chem.RWMol()

    # Add atoms from OBMol to RDKit Mol
    ob_atom_to_rd_atom = {}
    for ob_atom in openbabel.OBMolAtomIter(obmol):
        rd_idx = rdkit_mol.AddAtom(Chem.Atom(ob_atom.GetAtomicNum()))
        ob_atom_to_rd_atom[ob_atom.GetIdx()] = rd_idx

    # Add bonds from OBMol to RDKit Mol
    for ob_bond in openbabel.OBMolBondIter(obmol):
        rdkit_mol.AddBond(ob_atom_to_rd_atom[ob_bond.GetBeginAtomIdx()],
                          ob_atom_to_rd_atom[ob_bond.GetEndAtomIdx()],
                          Chem.BondType.values[ob_bond.GetBondOrder()])

    # Convert to a Chem.Mol object
    rdkit_mol = rdkit_mol.GetMol()

    # Set the 3D coordinates
    conformer = Chem.Conformer(rdkit_mol.GetNumAtoms())
    for ob_atom in openbabel.OBMolAtomIter(obmol):
        x, y, z = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()
        conformer.SetAtomPosition(ob_atom_to_rd_atom[ob_atom.GetIdx()], (x, y, z))

    rdkit_mol.AddConformer(conformer)

    return rdkit_mol


def set_radical_electrons(rd_mol, mol_charge):
    abs_mol_charge = abs(mol_charge)
    for idx, atom in enumerate(rd_mol.GetAtoms()):
        atom_num = atom.GetAtomicNum()
        if atom_num == 3:
            continue
        typical_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
        actual_valence = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
        if actual_valence < typical_valence:
            radical_charge_on_atom = int(typical_valence-actual_valence)
            if abs_mol_charge != 0:
                radical_charge_append = min(radical_charge_on_atom, abs_mol_charge)
                radical_charge_on_atom = radical_charge_on_atom - radical_charge_append
                abs_mol_charge = abs_mol_charge - radical_charge_append
            rd_mol.GetAtomWithIdx(idx).SetNumRadicalElectrons(radical_charge_on_atom)
            if atom_num == 6:
                rd_mol.GetAtomWithIdx(idx).SetProp('atomLabel', 'C')
    return rd_mol


def mol_pkl_2_mol_pics(pickle_path: str, output_dir: str = "molecule_images"):
    """
    Process a pickle file containing molecule information and generate PDF images for each molecule.

    Parameters:
    pickle_path (str): Path to the pickle file
    output_dir (str): Directory to save the output PDF files (default: "molecule_images")

    Returns:
    None
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load pickle file
    with open(pickle_path, 'rb') as f:
        f_data = pickle.load(f)
    smiles_fp_out= open(os.path.join(output_dir, "smiles.txt"), "w")
    # Process each molecule in the pickle file
    for idx, a_mol_info in enumerate(f_data):
        # Create XYZ file
        a_mol = a_mol_info.molecule
        a_new_atoms = Atoms(symbols=a_mol.labels, positions=a_mol.cart_coords)

        output_filename = os.path.join(output_dir, f"{idx}.pdf")
        if len(a_new_atoms)==1:
            plot_atom_to_pdf(atom_symbol=a_mol.labels[0], output_filename=output_filename, charge=a_mol.charge)
            continue

        xyz_filename = f"{idx}.xyz"
        a_new_atoms.write(xyz_filename)

        # Convert XYZ to OBMol
        obmol = xyz_2_db_mol(idx)

        # Convert OBMol to RDKit Mol
        rdkit_mol = obmol_to_rdkit_mol(obmol)

        # Set charge and spin multiplicity
        rdkit_mol.SetProp("_Name", f"Molecule {idx}")
        rdkit_mol.SetProp("Charge", str(a_mol.charge))
        rdkit_mol.SetProp("SpinMultiplicity", str(a_mol.spin_multiplicity))
        rdkit_mol = set_radical_electrons(rdkit_mol, a_mol.charge)

        # Plot molecule to PDF
        plot_molecule_to_pdf(rdkit_mol, output_filename, a_mol.charge)

        # Remove temporary XYZ file
        os.remove(xyz_filename)

        smiles= Chem.MolToSmiles(rdkit_mol)
        try:
            smiles_fp_out.write(f"{idx}\t{smiles}\t{a_mol.charge}\n")
        except:
            pass

    print(f"Processed {len(f_data)} molecules. Images saved in '{output_dir}' directory.")


def apply_species_filter(json_path: str):
    with open(json_path, 'r') as f:
        database_entries = json.load(fp=f)

    mol_entries = species_filter(
        database_entries,
        mol_entries_pickle_location='mol_entries.pickle',
        species_report='unfiltered_species_report.tex',
        species_decision_tree=no_species_decision_tree,
        coordimer_weight=lambda mol: (mol.penalty, mol.solvation_free_energy),
        generate_unfiltered_mol_pictures=False
    )

# Example usage
if __name__ == "__main__":
    mol_pkl_2_mol_pics(r"mol_entries.pickle", 'mol_pictures')