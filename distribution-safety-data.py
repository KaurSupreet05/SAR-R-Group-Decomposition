import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, spearmanr

file_path = "data-necrosis-ld20.dat"  # Replace with the actual file path
df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=[ "LD20_Necrosis"])


# Step 2: Inspect the data
print(df.head())
print(df.info())

# Step 2: Plot the distribution of LD20 cell viability
plt.figure(figsize=(10, 6))
sns.histplot(df['LD20_Necrosis'], kde=True, bins=30, color='skyblue', alpha=0.7)
plt.title('Distribution of LD20 Necrosis', fontsize=16)
plt.xlabel('LD20 Necrosis', fontsize=12)
plt.ylabel('Frequency', fontsize=12)
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()

# Step 3: Calculate correlations (if additional features are available)
# Example: Adding a synthetic feature for demonstration
df['Another_Feature'] = np.random.normal(loc=50, scale=15, size=500)
pearson_corr, _ = pearsonr(df['LD20_Necrosis'], df['Another_Feature'])
spearman_corr, _ = spearmanr(df['LD20_Necrosis'], df['Another_Feature'])

print(f"Pearson Correlation: {pearson_corr:.2f}")
print(f"Spearman Correlation: {spearman_corr:.2f}")

# Step 4: Visualize correlations
plt.figure(figsize=(8, 6))
sns.scatterplot(x='LD20_Necrosis', y='Another_Feature', data=df, alpha=0.6)
plt.title('Scatter Plot of LD20 Cell Viability vs Another Feature', fontsize=16)
plt.xlabel('LD20 Necrosis', fontsize=12)
plt.ylabel('Another Feature', fontsize=12)
plt.grid(axis='both', linestyle='--', alpha=0.7)
plt.show()

# Pairplot for pairwise relationships if there are more features
sns.pairplot(df, diag_kind='kde', corner=True)
plt.show()

