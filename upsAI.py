#!/bin/env python3

################################################################################

#upsAI
#Copyright 2025 by Elcid Aaron Pangilinan.
#All Rights Reserved

#This file is part of the upsAI distribution and governed by
#the MIT License. Please see the LICENSE file which should
#have been included in this package.

# Usage: python upsAI.py <input_fasta_path> <gene_location>

################################################################################ IMPORTS

import os
import sys
import pandas as pd
import numpy as np
from joblib import load
from Bio import SeqIO
from collections import Counter
from sklearn.svm import SVC, LinearSVC

################################################################################ PREDETERMINED VARIABLES

# Paths
output_dir = r"./results"
os.makedirs(output_dir, exist_ok=True)
cached_features_dir = r"./ref"

# Define allowed gene locations
valid_gene_locations = {"tag"}

# Best parameters for each gene location
best_params_dict = {
    'tag': {'C': 0.004497428727081967}
}

all_amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
all_tetrapeptides = [a1 + a2 + a3 + a4 for a1 in all_amino_acids for a2 in all_amino_acids for a3 in all_amino_acids for a4 in all_amino_acids]

################################################################################ FUNCTIONS

def load_query_sequences(query_fasta_path):
    """Loads sequences and their IDs from a FASTA file."""
    sequences = []
    for record in SeqIO.parse(query_fasta_path, "fasta"):
        sequences.append((record.id, str(record.seq)))
    return sequences

def translate_sequence(sequence):
    """Returns the uppercase version of the sequence."""
    return sequence.upper()

def calculate_tetrapeptide_counts(sequence):
    """Calculates the Tetrapeptide Composition (TC) of a given sequence."""
    sequence = translate_sequence(sequence)
    tetrapeptides = [sequence[i:i+4] for i in range(len(sequence) - 3)]
    tetrapeptide_count = Counter(tetrapeptides)

    # Create the tetrapeptide feature vector
    tetrapeptide_vector = [tetrapeptide_count[tetrapeptide] for tetrapeptide in all_tetrapeptides]
    return tetrapeptide_vector

def transform_to_features(sequences):
    """Transforms sequences into tetrapeptide count feature vectors."""
    feature_vectors = [calculate_tetrapeptide_counts(seq) for seq in sequences]
    return np.array(feature_vectors)

def run_prediction_on_query(gene_location, query_fasta_path):
    if gene_location not in valid_gene_locations:
        print(f"Error: Invalid gene location '{gene_location}'. Must be one of {valid_gene_locations}.")
        sys.exit(1)

    # Load precomputed training library
    cached_feature_path = f"{cached_features_dir}/{gene_location}.joblib"
    best_params = best_params_dict[gene_location]
    
    try:
        X_train, y_train = load(cached_feature_path)
    except FileNotFoundError:
        print(f"Error: Cached feature file not found for gene location '{gene_location}'.")
        sys.exit(1)

    # Train the SVM model with the best hyperparameters
    clf = LinearSVC(C=best_params["C"], dual='auto')
    clf.fit(X_train, y_train)

    # Load query sequences
    query_sequences = load_query_sequences(query_fasta_path)
    sequence_ids, sequences = zip(*query_sequences)

    # Transform query sequences into feature vectors
    X_query = transform_to_features(sequences)

    # Predict classes and generate output file name
    predicted_labels = clf.predict(X_query)
    fasta_basename = os.path.basename(query_fasta_path).rsplit('.', 1)[0]
    output_file = os.path.join(output_dir, f"{fasta_basename}_predictions.txt")
    
    # Save
    with open(output_file, "w") as f:
        f.write("Sequence_ID\tPredicted_Label\n")
        for seq_id, pred in zip(sequence_ids, predicted_labels):
            f.write(f"{seq_id}\t{pred}\n")
    print(f"Predictions saved to {output_file}")

################################################################################ MAIN

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python upsAI.py <input_fasta_path> <gene_location>")
        sys.exit(1)

    query_fasta_path = sys.argv[1]
    gene_location = sys.argv[2]

    if not os.path.exists(query_fasta_path):
        print(f"Error: File {query_fasta_path} does not exist.")
        sys.exit(1)

    run_prediction_on_query(gene_location, query_fasta_path)
