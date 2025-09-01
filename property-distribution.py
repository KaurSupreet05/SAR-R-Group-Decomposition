import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
import matplotlib.pyplot as plt

# File paths
input_file = "smiles_file.dat"  # Replace with your .dat file path
output_csv = "molecular_properties.csv"

# Read SMILES from .dat file
# Assuming each line in the file contains one SMILES string
with open(input_file, "r") as f:
    smiles_list = [line.strip() for line in f if line.strip()]

# Initialize property lists
molecular_weights = []
logPs = []
polar_surface_areas = []
rotatable_bonds = []
valid_smiles = []
invalid_smiles = []

# Compute properties
for smi in smiles_list:
    mol = Chem.MolFromSmiles(smi)
    if mol:  # Ensure valid molecule
        valid_smiles.append(smi)
        molecular_weights.append(Descriptors.MolWt(mol))
        logPs.append(Crippen.MolLogP(mol))
        polar_surface_areas.append(Descriptors.TPSA(mol))
        rotatable_bonds.append(Descriptors.NumRotatableBonds(mol))
    else:
        invalid_smiles.append(smi)

# Save invalid SMILES to a file (optional)
if invalid_smiles:
    with open("invalid_smiles.dat", "w") as f:
        f.writelines([f"{smi}\n" for smi in invalid_smiles])

# Create a DataFrame for valid SMILES and properties
data = pd.DataFrame({
    "SMILES": valid_smiles,
    "Molecular Weight (g/mol)": molecular_weights,
    "LogP": logPs,
    "Polar Surface Area (Å²)": polar_surface_areas,
    "Rotatable Bonds (count)": rotatable_bonds,
})

# Save DataFrame to CSV
data.to_csv(output_csv, index=False)
print(f"Molecular properties saved to {output_csv}")

# Plot histograms
properties = data.columns[1:]  # Exclude SMILES column
plt.figure(figsize=(14, 10))
for i, prop in enumerate(properties, 1):
    plt.subplot(2, 2, i)
    plt.hist(data[prop], bins=20, color='blue', edgecolor='black', alpha=0.7)
    plt.title(f"Distribution of {prop}")
    plt.xlabel(prop)
    plt.ylabel("Frequency")
    plt.grid(axis='y', alpha=0.75)

plt.tight_layout()
plt.show()

