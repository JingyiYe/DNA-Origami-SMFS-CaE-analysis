# AF-Cluster post-processing

This folder contains the Python script used for AF-Cluster post-processing and figure-source-data generation.

## File

- `AFcluster_postprocess.py`  
  Python script used to calculate mean pLDDT values, Cα RMSDs relative to a common reference structure, all-pair Cα–Cα distance features, principal-component coordinates and machine-readable output tables.

## Requirements

The required Python environment is provided in the root-level files:

```text
environment.yml
requirements.txt
```

To create the recommended conda environment from the repository root:

```bash
conda env create -f environment.yml
conda activate smfs-cae-analysis
```

## Input files

Expected inputs include AF-Cluster-generated PDB files, prediction metadata files when available, and a reference PDB structure for Cα RMSD calculation.

## Outputs

The script exports RMSD tables, mean pLDDT tables, pairwise-distance features, PCA coordinates, metadata tables and plots for downstream structural interpretation and figure preparation.

## Notes

Users should ensure that all compared structures contain consistent residue numbering and a common set of Cα atoms before RMSD or PCA analysis.

The downstream results should be interpreted together with experimental contour-length constraints, HMM-derived mechanical states and construct boundaries.

## Related files

The modified AF-Cluster notebook is located in:

```text
AF-Cluster_notebook/AFcluster.ipynb
```

The passive-mode HMM analysis scripts are located in:

```text
matlab/passive_hmm/
```
