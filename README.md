# upsAI: A Classifier for *Plasmodium falciparum* var Gene Upstream Groups

**upsAI** is a machine learning-based classifier that predicts the upstream group (upsA, upsB, upsC, or upsE) of *Plasmodium falciparum* var genes from tag sequences. This tool enables accurate and rapid assignment of var genes to upstream groups, an important step for understanding the parasite’s role in disease severity and immune evasion.

---

## About

The *PfEMP1* family of proteins, encoded by the hypervariable **var** gene family in *Plasmodium falciparum*, plays a central role in antigenic variation and virulence. Var genes can be grouped by their upstream regions (upsA, upsB, upsC, upsE), with known links between these groups and disease outcomes.
To address the need for an accurate classifier, we developed **upsAI** using a custom dataset of var gene sequences. Our model was developed from different k-mer decomposition and different machine learning algorithms.

---

### Usage

```bash
python upsAI.py <input_fasta_path> <gene_location>
```
