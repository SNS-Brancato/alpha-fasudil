import prolif as plf
import MDAnalysis as mda
from rdkit import Chem, DataStructs
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def load_universe(psf, dcd):
    """
    Loads the MD simulation data into a Universe object.
    
    Args:
        psf (str): Path to the PSF file.
        dcd (str): Path to the DCD trajectory file.
    
    Returns:
        MDAnalysis.Universe: Loaded universe object.
    """
    return mda.Universe(psf, dcd)

def generate_fingerprint(universe, ligand, protein, output_path):
    """
    Generates and saves the fingerprint of interactions between ligand and protein.
    
    Args:
        universe (MDAnalysis.Universe): The Universe object containing simulation data.
        ligand (MDAnalysis.Atoms): The ligand atoms in the simulation.
        protein (MDAnalysis.Atoms): The protein atoms in the simulation.
        output_path (str): Path to save the generated fingerprint.
    """
    fp = plf.Fingerprint()
    fp.run(universe.trajectory, ligand, protein)
    fp.to_pickle(output_path)
    return fp

def calculate_interaction_percentage(df):
    """
    Calculates the percentage of interactions for each interaction type.
    
    Args:
        df (pd.DataFrame): DataFrame containing interaction data.
    
    Returns:
        pd.DataFrame: Interaction percentages sorted in descending order.
    """
    return df.T.groupby('interaction').sum().T.astype(bool).mean().sort_values(ascending=False).to_frame(name="%").T * 100

def calculate_residue_interactions(df):
    """
    Calculates the percentage of residues most frequently interacting with the ligand.
    
    Args:
        df (pd.DataFrame): DataFrame containing interaction data.
    
    Returns:
        pd.DataFrame: Top 10 residues with the highest interaction percentage.
    """
    return (
        df.T.groupby(level=["ligand", "protein"])
        .sum()
        .T.astype(bool)
        .mean()
        .sort_values(ascending=False)
        .head(10)
        .to_frame("%")
        .T * 100
    )

def calculate_interaction_type(df, interaction_type):
    """
    Calculates the percentage of a specific interaction type.
    
    Args:
        df (pd.DataFrame): DataFrame containing interaction data.
        interaction_type (str): Interaction type to calculate (e.g., "PiStacking").
    
    Returns:
        pd.DataFrame: Interaction percentage for the specified interaction type.
    """
    return (
        df.xs(interaction_type, level="interaction", axis=1)
        .mean()
        .sort_values(ascending=False)
        .to_frame(name="%")
        .T * 100
    )

def calculate_tanimoto_similarity(fp):
    """
    Calculates the Tanimoto similarity matrix for the given fingerprint.
    
    Args:
        fp (plf.Fingerprint): The Fingerprint object.
    
    Returns:
        pd.DataFrame: Tanimoto similarity matrix.
    """
    bitvectors = fp.to_bitvectors()
    similarity_matrix = [
        DataStructs.BulkTanimotoSimilarity(bv, bitvectors) for bv in bitvectors
    ]
    return pd.DataFrame(similarity_matrix, index=fp.to_dataframe().index, columns=fp.to_dataframe().index)

def plot_heatmap(similarity_matrix):
    """
    Plots the heatmap for the Tanimoto similarity matrix.
    
    Args:
        similarity_matrix (pd.DataFrame): Tanimoto similarity matrix to plot.
    """
    fig, ax = plt.subplots(figsize=(3, 3), dpi=200)
    colormap = sns.diverging_palette(300, 145, s=90, l=80, sep=30, center="dark", as_cmap=True)
    sns.heatmap(
        similarity_matrix,
        ax=ax,
        square=True,
        cmap=colormap,
        vmin=0,
        vmax=1,
        center=0.5,
        xticklabels=5,
        yticklabels=5,
    )
    ax.invert_yaxis()
    plt.yticks(rotation="horizontal")
    fig.patch.set_facecolor("white")

def main(psf, dcd, output_fingerprint_path):
    """
    Main function to run the analysis pipeline.
    
    Args:
        psf (str): Path to the PSF file.
        dcd (str): Path to the DCD trajectory file.
        output_fingerprint_path (str): Path to save the fingerprint.
    """
    # Load simulation data
    u = load_universe(psf, dcd)
    
    # Select ligand and protein atoms
    lig = u.atoms.select_atoms('resname UNL')
    prot = u.atoms.select_atoms('protein')
    
    # Generate and save fingerprint
    fp = generate_fingerprint(u, lig, prot, output_fingerprint_path)
    
    # Convert fingerprint to DataFrame
    df = fp.to_dataframe()
    
    # Calculate interaction percentages
    interaction_percentage = calculate_interaction_percentage(df)
    
    # Calculate residue interactions
    residue_interactions = calculate_residue_interactions(df)
    
    # Calculate specific interaction types
    pi_stacking = calculate_interaction_type(df, "PiStacking")
    hydrophobic = calculate_interaction_type(df, "Hydrophobic")
    hba_acceptor = calculate_interaction_type(df, "HBAcceptor")
    hbd_donor = calculate_interaction_type(df, "HBDonor")
    vdW_contact = calculate_interaction_type(df, "VdWContact")
    cationic = calculate_interaction_type(df, "Cationic")
    cation_pi = calculate_interaction_type(df, "CationPi")
    
    # Calculate Tanimoto similarity
    tanimoto_similarity = calculate_tanimoto_similarity(fp)
    
    # Plot Tanimoto similarity heatmap
    plot_heatmap(tanimoto_similarity)
    
    # Return results for further analysis or display
    return {
        "interaction_percentage": interaction_percentage,
        "residue_interactions": residue_interactions,
        "pi_stacking": pi_stacking,
        "hydrophobic": hydrophobic,
        "hba_acceptor": hba_acceptor,
        "hbd_donor": hbd_donor,
        "vdw_contact": vdW_contact,
        "cationic": cationic,
        "cation_pi": cation_pi,
        "tanimoto_similarity": tanimoto_similarity
    }

if __name__ == "__main__":
    psf = '/home/taylan/Desktop/alpha/CNS_BBB/best_docked_10/md_prep/ready_to_protein/md_G610-0179pc.prmtop'
    dcd = '/home/taylan/Desktop/alpha/CNS_BBB/best_docked_10/md_prep/ready_to_protein/G610-0179_300ns.dcd'
    output_fingerprint_path = "/home/taylan/Desktop/alpha/CNS_BBB/best_docked_10/md_prep/ready_to_protein/G610-0179_fingerprint.pkl"
    
    results = main(psf, dcd, output_fingerprint_path)
    # You can print or save the results for further processing
    print(results)
