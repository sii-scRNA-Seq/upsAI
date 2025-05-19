# upsAI: A Classifier for *Plasmodium falciparum* var Gene Upstream Groups

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

| Option                | Description                                                         |
| --------------------- | ------------------------------------------------------------------- |
| `-i`, `--input-fasta` | Path to the input FASTA file containing var gene sequences          |
| `-m`, `--model-name`  | Prefix of the model file to use (e.g., `kmer_rf`)                   |
| `-o`, `--output-dir`  | Directory to save prediction results (default: `./results`)         |
| `-d`, `--model-dir`   | Directory containing trained models (default: `./models`)           |
| `--list-models`       | List all available models in the specified model directory and exit |
| `-v`, `--version`     | Show the current version of upsAI and exit                          |
```
