import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# 1. Load the CSV file
data_file = "LD20-bioactivity-smiles.csv"  # Replace with your file path
df = pd.read_csv(data_file)

# Ensure the DataFrame has columns: "SMILES" and "Bioactivity"
if not {"SMILES", "Bioactivity"}.issubset(df.columns):
    raise ValueError("CSV file must have 'SMILES' and 'Bioactivity' columns.")

# 2. Calculate 3D coordinates and shape descriptors
mols = []
shape_descriptors = []

for smiles in df["SMILES"]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        shape_descriptors.append([np.nan] * 3)
        mols.append(None)
        continue

    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        radius_of_gyration = rdMolDescriptors.CalcRadiusOfGyration(mol)
        inertial_shape_factor = rdMolDescriptors.CalcInertialShapeFactor(mol)
        
        # Using radius of gyration as a shape descriptor for volume
        shape_descriptors.append([radius_of_gyration, inertial_shape_factor])
        mols.append(mol)
    except Exception as e:
        print(f"Failed to process {smiles}: {e}")
        shape_descriptors.append([np.nan] * 2)
        mols.append(None)

# Add shape descriptors to the DataFrame
shape_df = pd.DataFrame(shape_descriptors, columns=["RadiusGyration", "InertialShape"])
df = pd.concat([df, shape_df], axis=1)

# Drop rows with NaN values
df = df.dropna()

# 3. Plot histograms for shape descriptors (Radius of Gyration and Inertial Shape Factor)
plt.figure(figsize=(14, 6))

# Radius of Gyration Histogram
plt.subplot(1, 2, 1)
plt.hist(df["RadiusGyration"], bins=30, color='blue', edgecolor='black', alpha=0.7)
plt.title("Distribution of Radius of Gyration")
plt.xlabel("Radius of Gyration")
plt.ylabel("Frequency")
plt.grid(True)

# Inertial Shape Factor Histogram
plt.subplot(1, 2, 2)
plt.hist(df["InertialShape"], bins=30, color='blue', edgecolor='black', alpha=0.7)
plt.title("Distribution of Inertial Shape Factor")
plt.xlabel("Inertial Shape Factor")
plt.ylabel("Frequency")
plt.grid(True)

plt.tight_layout()
plt.show()

