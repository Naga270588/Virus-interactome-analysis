
import pandas as pd
import networkx as nx
import numpy as np
from pathlib import Path
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Start the timer
start_time = time.time()

# Set the folder path
folder_path = Path("path_to_your_folder")

# Load the data
virus_proteins = pd.read_excel(folder_path / "VIRUS.xlsx") # Replace 'VIRUS.xlsx' with your file
ppi_data = pd.read_excel(folder_path / "Interactome_data.xlsx")

# Identify correct column names
C_v = virus_proteins['gene_id'].tolist()

# Build the graph from the PPI data
G = nx.from_pandas_edgelist(ppi_data, source=ppi_data.columns[0], target=ppi_data.columns[1])

# Filter out genes that are not in the PPI network
C_v = [gene for gene in C_v if gene in G.nodes]

# Define the compute_proximity function
def compute_proximity(G, C, T):
    min_dist_C = []
    min_dist_T = []
    
    for c in C:
        distances = [nx.shortest_path_length(G, c, t) for t in T if nx.has_path(G, c, t)]
        if distances:
            min_dist_C.append(min(distances))
    
    for t in T:
        distances = [nx.shortest_path_length(G, t, c) for c in C if nx.has_path(G, t, c)]
        if distances:
            min_dist_T.append(min(distances))
    
    sum_dist = sum(min_dist_C) + sum(min_dist_T)
    d_CT = sum_dist / (len(C) + len(T))
    
    return d_CT

# Define the permutation_test function
def permutation_test(G, C, T, num_permutations=1000):
    original_d_CT = compute_proximity(G, C, T)
    random_proximities = []
    
    all_genes = list(G.nodes)
    
    for _ in range(num_permutations):
        random_C = np.random.choice(all_genes, size=len(C), replace=False)
        random_T = np.random.choice(all_genes, size=len(T), replace=False)
        
        random_d_CT = compute_proximity(G, random_C, random_T)
        random_proximities.append(random_d_CT)
    
    # Compute mean and standard deviation of random proximities
    mean_d_r = np.mean(random_proximities)
    std_d_r = np.std(random_proximities)
    
    # Calculate Z-score for original proximity
    Z_d_CT = (original_d_CT - mean_d_r) / std_d_r
    
    # Calculate P-value
    p_value = np.mean([1 if random_prox < original_d_CT else 0 for random_prox in random_proximities])
    
    # Return all required values
    return Z_d_CT, p_value, original_d_CT, mean_d_r, std_d_r, random_proximities

# Define the function to process each file
def process_file(file_path):
    try:
        drug_target_data = pd.read_excel(file_path)
        
        combined_results = []
        
        for drug in drug_target_data['Drug'].unique():
            T_d = drug_target_data[drug_target_data['Drug'] == drug]['Target'].tolist()
            
            # Filter out genes that are not in the PPI network
            T_d = [gene for gene in T_d if gene in G.nodes]
            
            if not T_d or not C_v:
                print(f"No matching genes found in the PPI network for drug {drug}.")
                continue
            
            # Get the detailed results from the permutation test
            Z_d_CT, p_value, original_d_CT, mean_d_r, std_d_r, random_proximities = permutation_test(G, C_v, T_d)
            
            # Store combined results
            combined_results.append({
                'Drug': drug,
                'Z-score': Z_d_CT,
                'P-value': p_value,
                'Significant Proximity': Z_d_CT < -1.5 and p_value < 0.05,
                'Original d_CT': original_d_CT,
                'Mean d_r': mean_d_r,
                'Std d_r': std_d_r
            })
        
        return combined_results
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")
        return []

# File paths
file_paths = [
    folder_path / "Drug_Target_Part_1.xlsx"
]

# Initialize results storage
all_combined_results = []

# Use ThreadPoolExecutor to run the tasks in parallel
with ThreadPoolExecutor() as executor:
    future_to_file = {executor.submit(process_file, file_path): file_path for file_path in file_paths}
    
    for future in as_completed(future_to_file):
        file_path = future_to_file[future]
        try:
            combined_results = future.result()
            all_combined_results.extend(combined_results)
        except Exception as e:
            print(f"Error processing file {{file_path}}: {{e}}")

# Convert combined results to a DataFrame
combined_results_df = pd.DataFrame(all_combined_results)

# Save the combined results to an Excel file in the same folder
combined_results_df.to_excel(folder_path / "New1.xlsx", index=False)

# Calculate and print the runtime
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Analysis complete. Results saved to 'New1.xlsx'.")
print(f"Total runtime: {elapsed_time:.2f} seconds.")
