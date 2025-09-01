import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Example dataset: Replace with your actual file path
data = pd.read_csv('LD20-6-correlation-data.csv')  # Columns: ['Candidate', 'Cell_Viability', 'Necrosis', 'Cell_Death', 'Nuclei_Area']

# Compute the correlation matrix
correlation_matrix = data[['Cell_Viability', 'Necrosis', 'Cytoplasm_area', 'Cell_Death', 'Nuclei_Area', 'Bioactivity']].corr(method='pearson')

# Display the correlation matrix
print(correlation_matrix)

# Visualize the correlation matrix using a heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
plt.title('Correlation Matrix of LD20 Assays (in μM)')
plt.savefig('Correlation_Matrix_Heatmap-LD20-6.png', dpi=300, bbox_inches='tight')
plt.show()

