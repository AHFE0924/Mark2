# NDM-1 Mutation Hotspot Predictor

A Graph Neural Network (GNN) project that predicts future probable NDM-1 mutations using Protein Language Models and structural analysis.

## Overview

This project uses:
- ESM-2 embeddings (650M parameters) for evolutionary information
- ESMFold for structure prediction
- Ensemble Graph Neural Networks for mutation prediction
- Multi-modal scoring combining evolutionary and structural features

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Using the Jupyter Notebook
Open `notebook.ipynb` in Jupyter Lab or Jupyter Notebook:
```bash
jupyter notebook notebook.ipynb
```

### Using the Python Script
Run the standalone Python script:
```bash
python ndm_mutation_predictor.py
```

## Project Structure

- `notebook.ipynb` - Interactive Jupyter notebook with analysis
- `ndm_mutation_predictor.py` - Standalone Python script
- `requirements.txt` - Python dependencies

## Citation

This project implements methods from:
- ESM-2/ESMFold: Lin Z, et al. "Evolutionary-scale prediction of atomic-level protein structure with a language model." Science (2023).
- NDM-1 characterization: Yong D, et al. "Characterization of a new metallo-β-lactamase gene, blaNDM-1..." Antimicrob Agents Chemother (2009).
