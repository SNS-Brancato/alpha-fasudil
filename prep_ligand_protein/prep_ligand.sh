#!/bin/bash

# Array of molecule types
molecule_types=("compact" "s_shape" "elongated")
# Array of molecule numbers
molecule_numbers=(1 2 3)

# Loop over each molecule type and number
for molecule_type in "${molecule_types[@]}"; do
  for number in "${molecule_numbers[@]}"; do
    # Create the filenames based on the molecule type and number
    mol_name="${molecule_type}_${number}"
    mol2_file="${mol_name}_prep.mol2"
    frcmod_file="${mol_name}_prep.frcmod"
    lib_file="${mol_name}_prep.lib"
    prmtop_file="${mol_name}_prep.prmtop"
    rst7_file="${mol_name}_prep.rst7"
    tleap_input="${type}_${num}.tleap.in"

    # Run tleap for each molecule and number
    cat > "tleap_input" <<EOF
source leaprc.gaff2
M77 = loadmol2 $mol2_file
check M77
loadamberparams $frcmod_file
saveoff M77 $lib_file
saveamberparm M77 $prmtop_file $rst7_file
quit
EOF

    tleap -f "tleap_input"
    echo "Processed $mol_name"
  done
done
