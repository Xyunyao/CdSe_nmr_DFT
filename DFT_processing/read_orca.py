'''
scripts to read orca NMR parameters
input: file_path, nuclus (cap-sensitive)
'''

import numpy as np
import re
def extract_tensors_from_file(file_path, nucleus_filter=None):
    
    data_dict = {}
    
     # Define your regular expression patterns here
    nucleus_pattern = re.compile(r'Nucleus\s+(\d+\w+):')
    tensor_pattern = re.compile(r'Total shielding tensor \(ppm\):\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)', re.DOTALL)

    # Read the text from the provided file path
    with open(file_path, 'r') as file:
        txt = file.read()

    nucleus_matches = nucleus_pattern.finditer(txt)
    tensor_matches = tensor_pattern.finditer(txt)

    nuclei = []
    for match in nucleus_matches:
        nuclei.append(match.group(1))

    tensors = []
    for match in tensor_matches:
        tensors.append(match.groups())

    if nuclei and tensors:
        for nucleus, tensor_match in zip(nuclei, tensors):
            # Check if a nucleus filter is provided, and if so, only include matching nuclei
            if nucleus_filter is None or nucleus_filter in nucleus:
                tensor_values = [float(value) for value in tensor_match]
                tensor_array = np.array(tensor_values).reshape(3, 3)
                data_dict[nucleus] = tensor_array

        # Output data if matches are found
        if data_dict:
            print("Data extracted:")
            for nucleus, tensor_array in data_dict.items():
                print("Nucleus:", nucleus)
                print("Tensor array:")
                print(tensor_array)
                print()
        else:
            print("No matches found for the specified nuclei.")
    else:
        print("No matches found in the file.")

    return data_dict
