import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os

# Load CSV file into a pandas DataFrame
df = pd.read_csv("results/adjusted_alignments.csv")

# Ensure the 'results/alignments' directory exists
os.makedirs('results/alignments', exist_ok=True)

# Function to plot the gene locations on the chromosome and save the plot
def plot_and_save_gene_locations(df, gene_name):
    # Filter data for the selected gene
    gene_data = df[df['Gene'] == gene_name]
    
    if gene_data.empty:
        print(f"No data found for gene: {gene_name}")
        return
    
    # Extract the relevant information
    chromosome = gene_data['Chromosome'].iloc[0]  # Assuming the chromosome is the same for all rows
    min_start = gene_data['Start Location'].min()
    max_end = (gene_data['Start Location'] + gene_data['Length']).max()

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 2))

    # Draw the chromosome as a rectangle (extend +-1000 bp from min and max)
    ax.add_patch(patches.Rectangle(
        (min_start - 1000, 0), max_end - min_start + 2000, 1,
        linewidth=1, edgecolor='black', facecolor='lightgrey'))

    # Plot each genome's gene location as a colored block
    genomes = gene_data['Genome'].unique()
    colors = plt.cm.get_cmap('tab20', len(genomes))  # Use a color map for different genomes

    for i, genome in enumerate(genomes):
        genome_data = gene_data[gene_data['Genome'] == genome]
        for _, row in genome_data.iterrows():
            start = row['Start Location']
            end = start + row['Length']
            ax.add_patch(patches.Rectangle(
                (start - 1000, 0), end - start, 1,
                linewidth=1, edgecolor='black', facecolor=colors(i)))

            # Add labels for the genome
            ax.text(start + (end - start) / 2 - 1000, 0.5, genome, ha='center', va='center')

    # Set the axis limits
    ax.set_xlim(min_start - 1000, max_end + 1000)
    ax.set_ylim(-0.5, 1.5)

    # Add chromosome labels
    ax.text(min_start - 1000, 1.1, f"{chromosome}: {min_start - 1000} bp", ha='left', va='center')
    ax.text(max_end + 1000, 1.1, f"{max_end + 1000} bp", ha='right', va='center')

    # Remove axis ticks and labels
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Add title
    plt.title(f"Gene Location Visualization for {gene_name}")
    
    # Save the plot to a PDF
    pdf_path = f"results/alignments/{gene_name}_alignment.pdf"
    plt.savefig(pdf_path, format='pdf')
    plt.close()  # Close the plot to free up memory

    print(f"Saved plot for {gene_name} to {pdf_path}")

# Loop through all unique genes in the CSV and generate the plots
unique_genes = df['Gene'].unique()
for gene in unique_genes:
    plot_and_save_gene_locations(df, gene)
