import os
import pandas as pd
import networkx as nx
import numpy as np
from scipy.stats import ranksums
import time

# Start timing the analysis
start_time = time.time()

# Define the path to your data folder
data_folder = "path_to_your_folder"

# Load the interactome data (human PPI network)
interactome_file = os.path.join(data_folder, "Interactome_data.xlsx")
interactome_data = pd.read_excel(interactome_file)

# Construct the PPI network
G = nx.from_pandas_edgelist(interactome_data, 'gene_id_A', 'gene_id_B')

# Load the viral-interacting human protein set for VIRUS
denv2_file = os.path.join(data_folder, "VIRUS.xlsx") # Replace 'VIRUS.xlsx' with your file
denv2_data = pd.read_excel(virus_file)
denv2_proteins = set(virus_data['gene_id'].tolist())

# Print total number of unique VIRUS-interacting human proteins
print(f"Total VIRUS-interacting human proteins: {len(virus_proteins)}")

# Calculate centrality measures for the human PPI network
degree_dict = dict(G.degree())
closeness_dict = nx.closeness_centrality(G)
betweenness_dict = nx.betweenness_centrality(G)

# Create a DataFrame to hold the centrality values
centrality_df = pd.DataFrame({
    'protein': list(G.nodes),
    'degree': [degree_dict[p] for p in G.nodes],
    'closeness': [closeness_dict[p] for p in G.nodes],
    'betweenness': [betweenness_dict[p] for p in G.nodes],  # Original betweenness
    'Betweenness(*10^4)': [betweenness_dict[p] * 1e4 for p in G.nodes]  # Betweenness * 10^4
})

# Label proteins as VIRUS-interacting or not
centrality_df['virus_interacting'] = centrality_df['protein'].apply(lambda x: x in denv2_proteins)

# Apply log2 transformation to the degree and closeness centrality measures
centrality_df['log2_degree'] = np.log2(centrality_df['degree'] + 1)
centrality_df['log2_closeness'] = np.log2(centrality_df['closeness'] + 1)
centrality_df['log2_betweenness'] = np.log2(betweenness_df['betweenness'] + 1)

# Define output file paths
output_file = os.path.join(data_folder, "centrality_results_VIRUS.xlsx")

# Save the centrality DataFrame to an Excel file
centrality_df.to_excel(output_file, index=False)

print(f"Centrality results for VIRUS saved to {output_file}")

# Perform Wilcoxon rank-sum tests
degree_stat, degree_p = ranksums(
    centrality_df[centrality_df['virus_interacting']]['log2_degree'],
    centrality_df[~centrality_df['virus_interacting']]['log2_degree']
)

closeness_stat, closeness_p = ranksums(
    centrality_df[centrality_df['virus_interacting']]['log2_closeness'],
    centrality_df[~centrality_df['virus_interacting']]['log2_closeness']
)

betweenness_stat, betweenness_p = ranksums(
    centrality_df[centrality_df['virus_interacting']]['log2_betweenness'],
    centrality_df[~centrality_df['virus_interacting']]['log2_betweenness']
)

# Use 'Betweenness(*10^4)' column for betweenness test
betweenness_stat, betweenness_p = ranksums(
    centrality_df[centrality_df['virus_interacting']]['Betweenness(*10^4)'],
    centrality_df[~centrality_df['virus_interacting']]['Betweenness(*10^4)']
)

# Prepare the Wilcoxon test results
wilcoxon_results = pd.DataFrame({
    'Centrality': ['Degree', 'Closeness', 'Betweenness'],
    'Statistic': [degree_stat, closeness_stat, betweenness_stat],
    'p-value': [degree_p, closeness_p, betweenness_p]
})

# Save the Wilcoxon test results to an Excel file
wilcoxon_file = os.path.join(data_folder, "wilcoxon_test_results_VIRUS.xlsx")
wilcoxon_results.to_excel(wilcoxon_file, index=False)

print(f"Wilcoxon test results for VIRUS saved to {wilcoxon_file}")

# End timing the analysis and calculate the elapsed time
end_time = time.time()
elapsed_time = end_time - start_time

print(f"Total time taken for the analysis: {elapsed_time:.2f} seconds")
