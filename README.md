# Interaction of Alpha-Synuclein with Fasudil

This project aims to explore the interactions between the alpha-synuclein protein and the fasudil drug in the presence of Ca²⁺ ions.

## Project Structure

    docking/: Contains the necessary files to perform docking simulations using AutoDock Vina. These files include ligand and receptor structures, docking configurations, and scripts to generate docking results.

    prep_ligand_protein/: This directory holds files used for preparing the protein-drug complex for molecular dynamics (MD) simulations with Amber.

    Int_fasudil_synuclein.py: A Python script for analyzing the interactions between alpha-synuclein and fasudil, focusing on interaction fingerprints. 

    utils_interaction.py: A utility Python script that provides helper functions to analyze and process the interaction data.

## Workflow Overview

    Docking: Perform docking simulations to predict the binding mode of fasudil to alpha-synuclein using AutoDock Vina.

    Protein-Drug Complex Preparation: Prepare the protein-drug complex for MD simulations in Amber by generating the necessary topology and parameter files.

    Interaction Analysis: Use the Python scripts to analyze the interactions between the protein and the drug, including the identification of key interaction fingerprints.
