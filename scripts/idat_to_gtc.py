import argparse
import subprocess
import os

BASH_SOURCE = ". ~/.bash_profile"
CONDA_ACTIVATE = "conda activate genetics"

def init_workspace():
    
    os.system(BASH_SOURCE)
    os.environ["DOTNET_SYSTEM_GLOBALIZATION_INVARIANT"] = "1"
    
    return()

def convert_idats(input_bpm,
                  input_egt,
                  idat_dir,
                  output_dir,
                  gender_est_thresh,
                  num_threads):
    
    mycommand = ["iaap-cli",
                 "gencall",
                  input_bpm,
                  input_egt,
                  output_dir,
                 "--idat-folder",
                 idat_dir,
                 "--output-gtc",
                 "--gender-estimate-call-rate-threshold",
                 str(gender_est_thresh), # must convert float to str
                 "--num-threads",
                 str(num_threads)]       # must convert int to str
    
    subprocess.run(mycommand)
    
    return()

def main():
    # parse command line arguments.
    # require 2 input files: bpm (manifest), egt (cluster),
    # require 1 output file: gtc 
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-bpm",
                        type=str,
                        help="Path to the Illumina manifest file.")
    parser.add_argument("--input-egt",
                        type=str,
                        help="Path to the cluster file.")
    parser.add_argument("--idat-dir",
                        type=str,
                        help="Path to the directory where idat files are stored.")
    parser.add_argument("--num-threads",
                        type=int,
                        help="Number of parallel threads to run.")
    parser.add_argument("--gender-estimate-call-rate-threshold",
                        type=float,
                        default=-0.1,
                        help="Threshold for autosomal call rate for gender estimation. Off by default")
    parser.add_argument("--output-dir",
                        type=str,
                        help="Path to output directory containing the gtc files.")
    args=parser.parse_args()
    
    init_workspace()
    
    convert_idats(input_bpm=args.input_bpm,
                  input_egt=args.input_egt,
                  idat_dir=args.idat_dir,
                  output_dir=args.output_dir,
                  gender_est_thresh=args.gender_estimate_call_rate_threshold,
                  num_threads=args.num_threads)
    
    return()

if __name__ == "__main__":
    main()