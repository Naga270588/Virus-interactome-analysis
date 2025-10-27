import os
import itertools

# Function to execute the command for each combination
def execute_command(g1, g2,output_path):
    output_file = g1.split('.')[0] + "+" + g2.split('.')[0] + ".txt"
    output_file_path = os.path.join(output_path,output_file)
    command = "python separation.py --g1 {} --g2 {} -o {}".format(g1, g2,output_file_path)
    os.system(command)

# current working directory
cwd = os.getcwd()
os.makedirs('output')

# Get the list of text files in the folder
folder_path = cwd  # Update this with your folder path
text_files = [file for file in os.listdir(folder_path) if file.endswith(".txt")]

# Generate all combinations of file pairs
file_combinations = itertools.combinations(text_files, 2)

output_path = os.path.join(cwd,'output')
# Execute the command for each combination
for combination in file_combinations:

    g1, g2 = combination
    execute_command(g1, g2,output_path)
