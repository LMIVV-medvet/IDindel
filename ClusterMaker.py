import sys
import argparse
import os
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist
import matplotlib
import matplotlib.pyplot as plt
import csv

# Takes a presence/absence matrix, calculates the Jaccard distance index between each samples and clusters them. Outputs a dendrogram of the clustering, the matrix in a colored map and a csv file of the clusters.
# usage: python ClusterMaker.py [-h] --matrix <CSV file> --prefix <output prefix> [--cluster <CSV file>]
# Get arguments
parser = argparse.ArgumentParser(description="Takes a presence/absence matrix in csv format, calculates the Jaccard distance index between each samples and clusters them. Outputs a dendrogram of the clustering, the matrix in a colored map and a csv file of the clusters.")
parser.add_argument('-m','--matrix', metavar='csv', required=True, help="Path to the input matrix file in csv format. The csv is expected to have column and row identifiers")
parser.add_argument('-p','--prefix', metavar='string', required=True, help="Output prefix.")
parser.add_argument('-c','--cluster', metavar='csv', help="Optional, csv file to rename clusters containing a specific sequence. Format: Sequence,New_cluster_name")

parser.add_argument('-fs','--figsize', metavar='integer', type=int, default=100, help="Optional, figure size, default:100") 
parser.add_argument('-dw','--dlwidth', metavar='integer', type=int, default=4, help="Optional, dendrogram line width, default:4") 
parser.add_argument('-t','--threshold', metavar='float', type=float, default=0.5, help="Optional, distance threshold for clustering, default:0.5") 
parser.add_argument('-df','--dfont', metavar='integer', type=int, default=12, help="Optional, dendrogram font size, default:12")
parser.add_argument('-mf','--mfont', metavar='integer', type=int, default=8, help="Optional, matrix font size, default:8")

args = parser.parse_args()
matrix = args.matrix
prefix = args.prefix
DendrogramLineWidth = args.dlwidth
ClusteringThreshold = args.threshold
DendrogramFontSize = args.dfont
MatrixFontSize = args.mfont
MaxFigSize = args.figsize

RenameClusters = args.cluster

if not os.path.exists(matrix):
    sys.exit("Error: the input matrix file is not found.")
    
# Read the CSV file with headers for both columns and rows
df = pd.read_csv(matrix, index_col=0)

# Convert the DataFrame to a NumPy array
matrix = df.values

# Sanity check for the presence/absence matrix 
if not np.all((matrix == 0) | (matrix == 1)):
    sys.exit("Error: the input matrix file has values other than 1 and 0.")

# Calculate the distance matrix using Jaccard distance
distance_matrix = pdist(matrix, metric='jaccard')

# Perform hierarchical clustering using the average approach AKA UPGMA
linkage_matrix = linkage(distance_matrix, method='average')

# Define the threshold for clustering
clusters = fcluster(linkage_matrix, ClusteringThreshold, criterion='distance')

# Adjust figsize if the matrix is too big
num_rows = df.shape[0]
num_cols = df.shape[1]
if num_rows<MaxFigSize and num_cols<MaxFigSize:
    figsize = (num_rows, num_cols)
else:
    ToScale = max(num_rows,num_cols)
    ScalingFactor = MaxFigSize/ToScale
    figsize = (num_rows*ScalingFactor, num_cols*ScalingFactor)

# Create a figure with subplots for the dendrogram and the matrix on top of each other
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, gridspec_kw={'height_ratios': [1, 4]})

# Plot the dendrogram and get the order of the leaves
matplotlib.rcParams['lines.linewidth'] = DendrogramLineWidth
dendro = dendrogram(linkage_matrix, labels=df.index, ax=ax1, color_threshold=ClusteringThreshold)
ax1.set_title('')
ax1.set_xlabel('')
ax1.set_ylabel('Jaccard Distance', fontsize=DendrogramFontSize)
ax1.xaxis.set_ticks([])

# Reorder and transpose the matrix to fit the dendrogram
ordered_matrix = matrix[dendro['leaves'], :]
transposed_ordered_matrix = ordered_matrix.T

# Plot the transposed and reordered matrix
ax2.imshow(transposed_ordered_matrix, cmap='gray', aspect='auto')
ax2.set_xlabel('')
ax2.set_ylabel('Indels')
ax2.set_title('')
ax2.yaxis.set_ticks([])
ax2.set_xticks(np.arange(len(dendro['ivl'])))
ax2.set_xticklabels(dendro['ivl'], rotation=90, fontsize=MatrixFontSize)

# Save the figure to a file in SVG format
plt.tight_layout()
plt.savefig(prefix + '.svg', format='svg')
plt.close()

# Output clusters to a CSV file
cluster_output = pd.DataFrame({'Sequence': df.index, 'Cluster': clusters})
cluster_output.to_csv(prefix + '.csv', index=False)

print("Clustering complete\n")

##############################################################################
def update_clusters(input_csv, changes_csv):
    # Read the input CSV file
    with open(input_csv, mode='r') as infile:
        reader = csv.reader(infile)
        data = list(reader)

    # Read the changes CSV file
    with open(changes_csv, mode='r') as changesfile:
        changes_reader = csv.reader(changesfile)
        changes = list(changes_reader)

    # Create a dictionary to store the new cluster names for each sequence
    changes_dict = {}
    renamed_clusters = {}

    # Iterate over the changes and apply them to the input data
    for change in changes:
        sequence_to_find = change[0]
        new_cluster_name = change[1]

        # Find the cluster of the given sequence in the input data
        cluster_to_update = None
        for row in data:
            if row[0] == sequence_to_find:
                cluster_to_update = row[1]
                break

        # If the sequence is found, update the cluster name for all sequences in that cluster
        if cluster_to_update:
            if cluster_to_update in renamed_clusters.values():
                new_cluster_name = f"{new_cluster_name}_{renamed_clusters[cluster_to_update]}"
            renamed_clusters[cluster_to_update] = new_cluster_name
            renamed_clusters[new_cluster_name] = new_cluster_name
            
            for row in data:
                if row[1] == cluster_to_update or row[1] == renamed_clusters.get(cluster_to_update, cluster_to_update):
                    row[1] = new_cluster_name

    # Write the updated data back to the input CSV file
    with open(input_csv, mode='w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerows(data)

if RenameClusters:
    update_clusters(prefix + '.csv', RenameClusters)

