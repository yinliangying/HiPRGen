import os
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from ase.db import connect
from ase import Atoms
import openbabel

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

    # Process each molecule in the pickle file
    for idx, a_mol_info in enumerate(f_data):
        # Create XYZ file
        a_mol = a_mol_info.molecule
        a_new_atoms = Atoms(symbols=a_mol.labels, positions=a_mol.cart_coords)
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

        # Plot molecule to PDF
        output_filename = os.path.join(output_dir, f"{idx}.pdf")
        plot_molecule_to_pdf(rdkit_mol, output_filename, a_mol.charge)

        # Remove temporary XYZ file
        os.remove(xyz_filename)

    print(f"Processed {len(f_data)} molecules. Images saved in '{output_dir}' directory.")

# Example usage
if __name__ == "__main__":
    mol_pkl_2_mol_pics(r"mol_entries.pickle", 'molecule_images_v3')