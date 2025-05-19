#!/bin/env python3

################################################################################
# upsAI - A Classifier for Plasmodium falciparum var Gene Upstream Groups
# Copyright 2025 by Elcid Aaron Pangilinan.
# All Rights Reserved
################################################################################

import os
import sys
import time
import argparse
import numpy as np
from joblib import load
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
from sklearn.preprocessing import LabelEncoder

################################################################################ CONFIGURATION

# Tetrapeptides
all_amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
all_tetrapeptides = [a1 + a2 + a3 + a4 for a1 in all_amino_acids
                     for a2 in all_amino_acids
                     for a3 in all_amino_acids
                     for a4 in all_amino_acids]

################################################################################ FUNCTIONS

def load_query_sequences(query_fasta_path):
    """Loads sequences and IDs from a FASTA file."""
    sequences = []
    for record in SeqIO.parse(query_fasta_path, "fasta"):
        sequences.append((record.id, str(record.seq)))
    return sequences

def translate_sequence(sequence):
    """Translates nucleotide to amino acid if needed, otherwise returns uppercase amino acid sequence."""
    sequence = sequence.upper()

    # Heuristically determine if this is likely a nucleotide sequence
    nucleotide_bases = set("ACGTU")
    non_nucleotide_bases = set(sequence) - nucleotide_bases

    if len(non_nucleotide_bases) > 0:
        # Probably an amino acid sequence
        return sequence
    else:
        # Probably a nucleotide sequence: attempt translation
        print("Warning: Detected a nucleotide sequence. Translating to amino acids...")
        try:
            return str(Seq(sequence).translate(to_stop=True))
        except Exception as e:
            print(f"Error translating sequence: {e}")
            sys.exit(1)

def calculate_tetrapeptide_counts(sequence):
    """Calculates the Tetrapeptide Composition (TC) for a sequence."""
    sequence = translate_sequence(sequence)
    tetrapeptides = [sequence[i:i+4] for i in range(len(sequence) - 3)]
    tetrapeptide_count = Counter(tetrapeptides)
    tetrapeptide_vector = [tetrapeptide_count[tetrapeptide] for tetrapeptide in all_tetrapeptides]
    return tetrapeptide_vector

def transform_to_features(sequences):
    return np.array([calculate_tetrapeptide_counts(seq) for seq in sequences])

def run_prediction_on_query(model_name, query_fasta_path):
    """Detects model type and runs prediction using a trained model."""
    model_path = os.path.join(model_dir, f"{model_name}.joblib")
    if not os.path.exists(model_path):
        print(f"Error: No model found for '{model_name}'.\n"
              f"Expected model file: {model_path}")
        sys.exit(1)
    
    # Load model: could be (clf, label_encoder) tuple or just clf
    model_object = load(model_path)
    if isinstance(model_object, tuple):
        clf, label_encoder = model_object
        is_xgboost = True
    else:
        clf = model_object
        label_encoder = None
        is_xgboost = False

    # Load query data
    query_sequences = load_query_sequences(query_fasta_path)
    sequence_ids, sequences = zip(*query_sequences)
    X_query = transform_to_features(sequences)

    # Predict
    print(f"Running prediction using {'XGBoost' if is_xgboost else 'SVM'}...")
    start_time = time.time()
    predicted = clf.predict(X_query)
    end_time = time.time()

    # Decode if using XGBoost with LabelEncoder
    if is_xgboost and label_encoder:
        predicted = label_encoder.inverse_transform(predicted)

    print(f"Prediction completed in {end_time - start_time:.4f} seconds.")

    # Save predictions
    fasta_basename = os.path.basename(query_fasta_path).rsplit('.', 1)[0]
    output_file = os.path.join(output_dir, f"{fasta_basename}_{model_name}_predictions.txt")
    with open(output_file, "w") as f:
        f.write("Sequence_ID\tPredicted_Label\n")
        for seq_id, pred in zip(sequence_ids, predicted):
            f.write(f"{seq_id}\t{pred}\n")

    print(f"Predictions saved to {output_file}")

################################################################################ MAIN
if __name__ == "__main__":
    UPS_AI_VERSION = "1.0.0"

    parser = argparse.ArgumentParser(
        description="upsAI - Classifier for Plasmodium falciparum var Gene Upstream Groups"
    )

    # Arguments
    parser.add_argument(
        "-i", "--input-fasta", 
        required=False,
        help="Path to input FASTA file"
    )
    parser.add_argument(
        "-m", "--model-name", 
        required=False,
        help="Gene model name (prefix of model file)"
    )

    parser.add_argument(
        "-o", "--output-dir", 
        default="./results", 
        help="Directory to save prediction results (default: ./results)"
    )
    parser.add_argument(
        "-d", "--model-dir", 
        default="./models", 
        help="Directory containing trained models (default: ./models)"
    )
    parser.add_argument(
        "--list-models",
        action="store_true",
        help="List available models in the specified model directory and exit"
    )
    parser.add_argument(
        "--version", "-v", 
        action="version",
        version=f"upsAI version {UPS_AI_VERSION}"
    )

    # Parse the arguments
    args = parser.parse_args()

    # Handle --list-models option
    if args.list_models:
        if not os.path.exists(args.model_dir):
            print(f"Error: Model directory '{args.model_dir}' does not exist.")
            sys.exit(1)

        model_files = sorted(f for f in os.listdir(args.model_dir) if f.endswith(".joblib"))
        if not model_files:
            print(f"No models found in directory: {args.model_dir}")
        else:
            print("Available Models:")
            for f in model_files:
                model_name = f.replace(".joblib", "")
                print(f" - {model_name}")
        sys.exit(0)

    # If --list-models was not used, input fasta and model name must be specified
    if not args.input_fasta or not args.model_name:
        print("Error: Input FASTA and model name must be specified unless --list-models is used.")
        sys.exit(1)

    # Set directories
    output_dir = args.output_dir
    model_dir = args.model_dir

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Check input file exists
    if not os.path.exists(args.input_fasta):
        print(f"Error: Input FASTA file '{args.input_fasta}' does not exist.")
        sys.exit(1)

    # Run prediction
    run_prediction_on_query(args.model_name, args.input_fasta)