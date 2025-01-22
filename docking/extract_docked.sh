#!/bin/bash

# Input PDBQT file
input_file="compact.pdbqt"

# Output PDB file prefix
output_prefix="compact_pose_"

# Convert the PDBQT file into multiple PDB files
/home/taylan/openbabel-openbabel-2-4-0/build/bin/obabel -i pdbqt "$input_file" -o pdb -O "${output_prefix}%d.pdb" --separate


