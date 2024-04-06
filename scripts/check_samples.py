import os
import re
import argparse

# Find all idat files from the input directories and subdirectories
def find_idat_files(directory):
    idat_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".idat"):
                idat_files.append(file)
    return idat_files

# Create a list of all sample ids
def sample_list(idat_files):
    samples=set()
    for idat_file in idat_files:
        sample = re.sub("_(Red|Grn)\\.idat", "", idat_file)
        samples.add(sample)
    return samples
            
# Check if each sample has bot Red and Grn .idat files                
def check_files(samples, idat_files):
    sample_validity = {}
    sample_validity["samples"] = []
    sample_validity["incomplete_samples"] = []
    for sample in samples:
        red_file = f"{sample}_Red.idat"
        grn_file = f"{sample}_Grn.idat"
        if red_file in idat_files and grn_file in idat_files:
            sample_validity["samples"].append(sample)
        else:
            sample_validity["incomplete_samples"].append(sample)
    return sample_validity

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--idat-dir",
                        type=str,
                        help="Path to directory with IDAT files.")
    parser.add_argument("--output-dir",
                        type=str, 
                        help="Filepath for output txt files.")
    args = parser.parse_args()
        
    input_directory = args.idat_dir
    idat_files = find_idat_files(input_directory)
    samples = sample_list(idat_files)
    sample_validity = check_files(samples, idat_files)
    valid_samples = sample_validity["samples"]
    not_valid_samples = sample_validity["incomplete_samples"]

    with open(f"{args.output_dir}/samples.txt", "w") as output_file:
        for sample in valid_samples:
            output_file.write(f"{sample}\n")
    
    with open(f"{args.output_dir}/incomplete_samples.txt", "w") as output_file:
        for sample in not_valid_samples:
            output_file.write(f"{sample}\n")
        
if __name__ == "__main__":
    main()


