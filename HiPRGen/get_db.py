import json
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import cclib
from ase.units import Hartree
from tqdm import tqdm
from ase.atoms import Atoms
from ase.db import connect


def gau_readout(filepath):
    info = cclib.io.ccread(filepath)
    a_new_atoms = Atoms(symbols=info.atomnos, positions=info.atomcoords[0])
    data = {'scfenergies': info.scfenergies[0], 'mult': info.mult, 'charge': info.charge,
            'zpve': info.zpve*Hartree, 'atomcharges': info.atomcharges,
            'enthalpy': (info.enthalpy * Hartree - info.scfenergies)[0], 'entropy': info.entropy * Hartree}
    return info.scfenergies, info.entropy, info.enthalpy, data, a_new_atoms


with open(r'libe.json') as f:
    f_data = json.load(f)

gau_src_path = r'/root/parse_gau_llibe/0708_libe'

# with open(r'sample.json') as f:
#     f_data = json.load(f)
#
# gau_src_path = r'F:\new_develop\test\HiPRGen\test\0703\test_nbo\sample'

cwd_ = os.getcwd()
electronic_energy_diff = []
enthalpy_diff = []
entropy_diff = []
failed_count = 0

with connect('dump.db') as dump_db:
    for a_mol_info in tqdm(f_data, desc="Processing molecules", unit="mol"):
        an_id = a_mol_info['molecule_id']
        os.chdir(gau_src_path)
        os.chdir(str(an_id))
        old_ee = a_mol_info['thermo']['eV']['electronic_energy']
        old_entropy = a_mol_info['thermo']['eV']['total_entropy']
        old_enthalpy = a_mol_info['thermo']['eV']['total_enthalpy']
        try:
            new_scfenergies, new_entropy, new_enthalpy, data, a_new_atoms = gau_readout('gau.log')
            electronic_energy_diff.append(new_scfenergies - old_ee)
            entropy_diff.append(new_entropy * Hartree - old_entropy)
            enthalpy_diff.append(new_enthalpy * Hartree - new_scfenergies - old_enthalpy)
            data['molecule_id'] = an_id
            dump_db.write(a_new_atoms, data=data)

        except Exception as e:
            print(f"Failed to process molecule ID {an_id}: {e}")
            failed_count += 1

os.chdir(cwd_)

# Convert lists to numpy arrays for easier manipulation
electronic_energy_diff = np.array(electronic_energy_diff).squeeze(-1)
entropy_diff = np.array(entropy_diff)
enthalpy_diff = np.array(enthalpy_diff).squeeze(-1)

# Visualization
sns.set(style="whitegrid", rc={"axes.grid": False})

plt.figure(figsize=(12, 8))

# Plot electronic energy difference
plt.subplot(3, 1, 1)
sns.histplot(electronic_energy_diff, kde=True, color=sns.color_palette("muted")[2], legend='Electronic')
plt.xlabel('Energy Difference (eV)')
plt.ylabel('Counts')

# Plot entropy difference
plt.subplot(3, 1, 2)
sns.histplot(entropy_diff, kde=True, color=sns.color_palette("pastel")[0], legend='Entropy')
plt.xlabel('Entropy Difference (eV)')
plt.ylabel('Counts')

# Plot enthalpy difference
plt.subplot(3, 1, 3)
sns.histplot(enthalpy_diff, kde=True, color=sns.color_palette("bright")[4], legend='Enthalpy')
plt.xlabel('Enthalpy Difference (eV)')
plt.ylabel('Counts')

plt.tight_layout()
plt.savefig('e_compare.png', dpi=600)


print("All plots have been saved.")
# Store the results in a JSON file
results = {
    "electronic_energy_diff": electronic_energy_diff.tolist(),
    "entropy_diff": entropy_diff.tolist(),
    "enthalpy_diff": enthalpy_diff.tolist()
}

with open('results.json', 'w') as f:
    json.dump(results, f)

print(f"Number of failed reads: {failed_count}")