# upsAI: A high-accuracy machine learning classifier for predicting *Plasmodium falciparum* var gene upstream groups 

**upsAI** is a machine learning-based classifier that predicts the upstream group (upsA, upsB, upsC, or upsE) or localization (upsA, upsB subtelomeric, upsB/upsC internal, or upsE) of *Plasmodium falciparum* var genes from tag, cassette, exon 1, or entire var gene sequences. This tool enables accurate and rapid assignment of var genes to upstream groups, an important step for understanding the parasite’s role in disease severity and immune evasion.

---

## About

The *PfEMP1* family of proteins, encoded by the hypervariable **var** gene family in *Plasmodium falciparum*, plays a central role in antigenic variation and virulence. Var genes can be grouped by their upstream regions (upsA, upsB, upsC, upsE), with known links between these groups and disease outcomes.
To facilitate accurate and scalable categorization of var genes, **upsAI** uses machine learning models trained on k-mer representations of upstream sequences.

---

## Installation

Clone this repository and install the dependencies:

```bash
git clone https://github.com/sii-scRNA-Seq/upsAI.git
cd upsAI
pip install -r requirements.txt
tar -xvzf ./models/tag_abc_linear.tar.gz ./models/tag_intsub_linear.tar.gz -C ./models/
```

---

## Usage

### Basic Command

```bash
python upsAI.py -i <input_fasta_path> -m <model_name>
```

### Full Options

```bash
python upsAI.py [OPTIONS]
```

| Option                | Description                                                         |
| --------------------- | ------------------------------------------------------------------- |
| `-i`, `--input-fasta` | Path to the input FASTA file containing var gene sequences          |
| `-m`, `--model-name`  | Prefix of the model file to use (e.g., `kmer_rf`)                   |
| `-o`, `--output-dir`  | Directory to save prediction results (default: `./results`)         |
| `-d`, `--model-dir`   | Directory containing trained models (default: `./models`)           |
| `--list-models`       | List all available models in the specified model directory and exit |
| `-v`, `--version`     | Show the current version of upsAI and exit                          |


**Note:**
- Only a subset of pre-trained models is included in this GitHub repository.
- Additional models can be downloaded from 10.5281/zenodo.15462399
- Downloaded models are provided as .tar.gz archives and must be decompressed before use.

Example:
```bash
tar -xvzf tag_abc_linear.tar.gz -C ./models/
```

---

## Model Performance

The table below summarizes the classification accuracy of the top-performing models for different input regions and classification tasks. Models were evaluated based on prediction accuracy of **upstream group (ups type)** and **subcellular localization** using various feature regions (e.g., full var gene, exon 1, cassette, tag).

### Upstream Group (ups type) Classification (abc/e: upsA, upsB, upsC, upsE)

| Input Region | Accuracy | Top Models                                       |
|--------------|----------|--------------------------------------------------|
| tag          | 0.85     | SVM Poly, SVM RBF                                |
| cassette     | 0.87     | SVM Linear                                       |
| exon 1       | 0.90     | SVM Sigmoid, XGBoost                             |
| var          | 0.92     | SVM Linear, SVM RBF, SVM Sigmoid                 |

### Subcellular Localization Prediction (intsub: upsA, upsB subtelomeric, upsB/upsC internal, upsE)

| Input Region | Accuracy | Top Models                                          |
|--------------|----------|-----------------------------------------------------|
| tag          | 0.78     | SVM Poly, SVM RBF, Random Forest                    |
| cassette     | 0.82     | SVM RBF                                             |
| exon 1       | 0.84     | SVM Linear, SVM RBF, XGBoost                        |
| var          | 0.92     | SVM Linear, SVM Poly, SVM RBF, SVM Sigmoid, XGBoost |

**Note:** var-based models generally offer the highest accuracy across both classification tasks.