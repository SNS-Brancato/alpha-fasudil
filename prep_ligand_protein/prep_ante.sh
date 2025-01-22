# Loop over each molecule type and number
for molecule_type in "${molecule_types[@]}"; do
  for number in "${molecule_numbers[@]}"; do
    # Create the filenames based on the molecule type and number
    mol_name="${molecule_type}_${number}"
    pdb_file="${mol_name}_prep.pdb"
    mol2_file="${mol_name}_prep.mol2"
    frcmod_file="${mol_name}_prep.frcmod"

    # Run antechamber to generate the mol2 file
    antechamber -i "$pdb_file" -fi pdb -o "$mol2_file" -fo mol2 -nc 1 -c bcc -at gaff2

    # Run parmchk2 to generate the frcmod file
    parmchk2 -i "$mol2_file" -f mol2 -o "$frcmod_file"

    echo "Processed $mol_name"
  done
done
