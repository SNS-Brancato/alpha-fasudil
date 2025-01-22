import MDAnalysis as mda
import prolif as plf
import matplotlib.pyplot as plt

def generate_fingerprint(psf_path, dcd_path, pickle_path, plot_path):
    """
    Generate a fingerprint for a ligand-protein interaction and save the fingerprint data and plot.

    Parameters:
    - psf_path (str): Path to the PSF file.
    - dcd_path (str): Path to the DCD file.
    - pickle_path (str): Path to save the fingerprint as a pickle file.
    - plot_path (str): Path to save the fingerprint barcode plot.
    """
    # Load the universe
    u = mda.Universe(psf_path, dcd_path)

    # Select ligand and protein atoms
    lig = u.atoms.select_atoms('resname M77')
    prot = u.atoms.select_atoms('protein')

    # Run the fingerprint calculation
    fp = plf.Fingerprint()
    fp.run(u.trajectory, lig, prot)

    # Save fingerprint data to a pickle file
    fp.to_pickle(pickle_path)

    # Convert fingerprint data to a DataFrame and transpose
    df = fp.to_dataframe().T

    # Plot and save the barcode image
    fp.plot_barcode()
    plt.savefig(plot_path)
    plt.close()  # Close the plot to free memory

# Example usage:
generate_fingerprint(
    psf_path="/mnt/hdd/fasudil/md_compact_3.prmtop",
    dcd_path="/mnt/hdd/fasudil/compact_3_500ns.dcd",
 pickle_path="/home/taylan/Desktop/alpha/docking/fasudil/prep_fasudil/md_compact_3_500ns.pkl",
   plot_path="/home/taylan/Desktop/alpha/docking/fasudil/prep_fasudil/md_compact_3_500ns.png"
)
