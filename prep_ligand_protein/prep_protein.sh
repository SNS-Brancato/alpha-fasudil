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
    tleap_input="${mol_name}.tleap.in"

    # Generate tleap input script
    cat > "$tleap_input" <<EOF
source leaprc.protein.ff99SBdisp
source leaprc.water.tip4pd_disp
source leaprc.gaff2
loadamberparams frcmod.taylan_ion
loadamberparams $frcmod_file
loadoff $lib_file

# Load the protein and drug structures
alphasyn = loadPdb "${molecule_type}.pdb"
drug = loadmol2 "$mol2_file"
complex = combine {alphasyn drug}

# Solvate the complex and add ions
solvateOct complex TIP4PDBOX 14.0
addIons complex CA 50
addIons complex Cl- 90

# Save output files
saveAmberParm complex m150/md_${mol_name}.prmtop m150/md_${mol_name}.inpcrd
savepdb complex m150/md_${mol_name}.pdb
quit
EOF

    # Run tleap for each molecule
    tleap -f "$tleap_input"
    echo "Processed $mol_name"
  done
done
