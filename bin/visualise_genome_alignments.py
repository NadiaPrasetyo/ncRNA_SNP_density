import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import numpy as np

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

    # Ensure the chromosome rectangle extends 1000 bp on both sides
    chromosome_start = min_start - 1000
    chromosome_end = max_end + 1000

    # Define distinct colors for each genome (more variation)
    genomes = gene_data['Genome'].unique()
    num_genomes = len(genomes)
    colors = plt.cm.tab20c(np.linspace(0, 1, num_genomes))  # Use 'tab20c' for a broader color palette

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 2))

    # Draw the chromosome as a rectangle (extend +-1000 bp from min and max)
    ax.add_patch(patches.Rectangle(
        (chromosome_start, 0.3), chromosome_end - chromosome_start, 0.3,  # Adjusted chromosome rectangle
        linewidth=1, edgecolor='black', facecolor='lightgrey'))

    # Plot each genome's gene location as a colored block
    for i, genome in enumerate(genomes):
        genome_data = gene_data[gene_data['Genome'] == genome]
        for _, row in genome_data.iterrows():
            start = row['Start Location']
            end = start + row['Length']
            ax.add_patch(patches.Rectangle(
                (start - 1000, 0.3), end - start, 0.3,  # Gene block moved upwards to y=0.1
                linewidth=1, edgecolor=colors[i], facecolor=colors[i]))

            # Add labels for the genome with 90-degree rotation and slight downward offset
            ax.text(start + (end - start) / 2 - 1000, -0.2, genome, ha='center', va='center', rotation=90)

    # Set the axis limits
    ax.set_xlim(chromosome_start, chromosome_end)
    ax.set_ylim(-0.5, 1.5)

    # Add chromosome labels
    ax.text(chromosome_start, 1.1, f"{chromosome}: {chromosome_start} bp", ha='left', va='center')
    ax.text(chromosome_end, 1.1, f"{chromosome_end} bp", ha='right', va='center')

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
