🧬 Network-Based Characterization of Human–Virus Interactions and Drug Repurposing

📖 Overview

This repository presents a network biology framework to explore human–virus protein–protein interactions (PPIs), elucidate viral disease mechanisms, and identify potential antiviral drugs and synergistic drug combinations.
It integrates all major network-based analyses — virus–human PPI centrality, network separation, target overlap, and drug–target proximity — into a single workflow.

🧩 Key Analyses
1. Virus–Human PPI Centrality

Builds a human PPI network and labels virus-interacting proteins.

Calculates Degree, Closeness, and Betweenness centrality.

Compares virus-interacting vs non-interacting proteins using Wilcoxon rank-sum tests.

2. Network Separation

Quantifies the topological distance between two sets of genes or proteins (virus or drug targets).

Negative separation → overlapping modules (shared pathways).

Positive separation → distinct functional modules.

Automated scripts calculate pairwise separation across multiple sets.

3. Virus Target Overlap

Measures the similarity between virus target sets using:

Jaccard Index

Overlap Coefficient

Identifies shared and unique host dependencies among different viruses.

4. Drug Repurposing via Network Proximity

Computes network proximity between virus-interacting proteins and drug target proteins.

Uses degree-preserving permutation tests to calculate Z-scores and p-values.

Prioritizes drugs with Z < −1.5 and p < 0.05 as potential antivirals.

Predicts synergistic drug pairs with negative separation (s < 0).

🧪 Highlights

Integrates 222,433 experimentally supported human–virus PPIs from 10 databases (HVIDB, HPIDB, VirusHostNet, PHISTO, NCBI, DenHunt, VirusMint, MINT, IntAct, BioGRID).

Analyzes 20 medically important viruses.

Identifies 19 prioritized repurposable drugs and 9 potential synergistic combinations.

Validated by drug–gene signature and transcriptomics data.

⚙️ Requirements

Python 3.8+

Dependencies: pandas, numpy, networkx, scipy, openpyxl

Install dependencies:

pip install pandas numpy networkx scipy openpyxl

▶️ Usage

Prepare Input Files

Interactome_data.xlsx → Human PPI network

<VIRUS>.xlsx → Virus-interacting human proteins

Drug_Target_Part_*.xlsx → Drug–target interactions

Run Analyses

```
python centrality_analysis.py
python separation.py
python Drug_Proximity_Analysis.py
```


Each script outputs results as .xlsx or .csv files for further visualization.

📊 Outputs

Centrality measures and statistical comparisons

Network separation and overlap scores

Virus–virus target similarity metrics

Drug proximity scores with Z-scores and p-values

Ranked antiviral drug and drug-combination candidates
