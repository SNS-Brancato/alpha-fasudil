# Interaction of Alpha-Synuclein with Fasudil

This project aims to explore the interactions between the alpha-synuclein protein and the fasudil drug in the presence of Ca²⁺ ions.

## Environment Setup

The environment.yml file contains the necessary packages to reproduce the analysis. Use the following commands to set up the environment:
```
conda env create -f environment.yml
conda activate alpha-fasudil-env
```


## Project Structure

docking/: Contains the necessary files to perform docking simulations using AutoDock Vina.

prep_ligand_protein/: This directory holds files used for preparing the protein-drug complex for molecular dynamics simulations with Amber.

Int_fasudil_synuclein.py: A Python script for analyzing the interactions between alpha-synuclein and fasudil, focusing on interaction fingerprints. 

utils_interaction.py: A utility Python script that provides helper functions to analyze and process the interaction data.

## Workflow Overview

Docking: Perform docking simulations to predict the binding mode of fasudil to alpha-synuclein using AutoDock Vina.

Protein-Drug Complex Preparation: Prepare the protein-drug complex for MD simulations in Amber by generating the necessary topology and parameter files.

Interaction Analysis: Use the Python scripts to analyze the interactions between the protein and the drug, including the identification of key interaction fingerprints.


![Image](https://github.com/user-attachments/assets/2e397ae9-07da-4a05-89ca-b3a524e8beb3)



