
import numpy as np
import argparse
import os
import csv

# Function to read coordinates from a file and extract Cd2+ and Se2- coordinates
def read_coordinates(file_path):
    cd_coords = []
    se_coords = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) < 4:
                continue  # Skip lines that don't have enough data
            try:
                atom_type = parts[0]
                coords = list(map(float, parts[1:4]))  # Attempt to convert to floats
                if atom_type == 'Cd':
                    cd_coords.append(coords)
                elif atom_type == 'Se':
                    se_coords.append(coords)
            except ValueError:
                # Skip lines that cannot be converted (like headers or comments)
                continue
    return cd_coords, se_coords

# Function to calculate the Euclidean distance between two points in 3D space
def euclidean_distance(coord1, coord2):
    return np.sqrt(np.sum((np.array(coord1) - np.array(coord2)) ** 2))

# Function to identify the position of a Cd2+ ion based on nearby Se2- ions
def identify_position(cd_coord, se_coords, threshold):
    connections = []
    for se_coord in se_coords:
        if euclidean_distance(cd_coord, se_coord) <= threshold:
            connections.append(tuple(se_coord))  # Store Se2- coordinates as tuples
    num_connections = len(connections)
    
    if num_connections == 1:
        return 'Vertex', connections
    elif num_connections == 2:
        return 'Edge', connections
    elif num_connections == 3:
        return 'Face', connections
    elif num_connections == 4:
        return 'Center', connections
    else:
        return 'Unknown', connections

# Function to write Cd2+ coordinates and position info to a CSV file
def write_cd_positions_with_info_csv(cd_coords, se_coords, output_file, threshold):
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write the header
        csvwriter.writerow(['Cd_X', 'Cd_Y', 'Cd_Z', 'Position', 'Se Connections'])
        
        for cd_coord in cd_coords:
            position, connections = identify_position(cd_coord, se_coords, threshold)
            # Format Se connections as strings
            se_connection_str = '; '.join([f"({conn[0]}, {conn[1]}, {conn[2]})" for conn in connections])
            csvwriter.writerow([cd_coord[0], cd_coord[1], cd_coord[2], position, se_connection_str])

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Extract Cd positions and their connectivity information from an XYZ file.")
    parser.add_argument("input_file", help="Path to the input XYZ file.")
    parser.add_argument("--cdse_threshold", type=float, default=2.8, help="Distance threshold for Cd-Se connectivity (default: 2.8).")
    args = parser.parse_args()

    # Get the input file path and extract the directory
    input_file = args.input_file
    input_dir = os.path.dirname(input_file)

    # Read Cd2+ and Se2- coordinates from the file
    cd_coords, se_coords = read_coordinates(input_file)

    # Define the output CSV file path in the same folder as input, with a modified name
    output_file = os.path.join(input_dir, "cd_positions_with_info.csv")

    # Write Cd2+ positions and position info to the CSV file
    write_cd_positions_with_info_csv(cd_coords, se_coords, output_file, args.cdse_threshold)

    print(f"Cd2+ positions and connectivity information written to {output_file}")

if __name__ == "__main__":
    main()


