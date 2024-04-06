import os
import re

# Find all idat files from the input directories and subdirectories
def find_idat_files(directory):
    idat_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".idat"):
                idat_files.append(file)
    return idat_files