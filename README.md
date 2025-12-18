# Do Regions with More Brokers Produce More Atypical Technology Combinations?
## Investigating the Role of Inventor Brokerage and Regional Innovation

### Course: Economic Complexity – Social Data Science (MSc)
**Institution:** Corvinus University of Budapest  
**Year:** 2025

### Team Members:
-   Patrik Bruinsma
-   Ákos Virág

### Project Overview
This repository contains the research pipeline for a project analyzing the relationship between network structure and regional innovation novelty. Innovation often emerges from the recombination of existing technologies, yet only a small fraction of inventions relies on atypical (rare) combinations of knowledge domains.

The primary objective of this research is to connect individual-level brokerage mechanisms with regional-level outcomes. We investigate whether regions with a higher density of "inventor brokers" (individuals bridging structural holes or occupying central positions) in their co-inventor networks generate a higher share of atypical technology combinations. The study employs a **Two-Way Fixed Effects (TWFE)** model across European NUTS 2 regions (1979–2023) to isolate the effect of brokerage from regional heterogeneity and temporal shocks.

### Key Findings
1.  **Short-Term vs. Long-Term:** A higher share of brokers does not drive immediate "spikes" in atypical patenting (flow), but significant results are found for the accumulated **stock** of atypical capabilities.
2.  **Network of Places:** We utilize the "Network of Places" transformation to correct for spurious density in large inventor teams.
3.  **The "Smart Spec" Interaction:** We find a significant positive interaction between brokerage (Mean Betweenness) and **Relatedness Density**, suggesting brokers are most effective in regions with coherent, related knowledge bases.

### Repository Structure
```text
.
├── data/
│   └── Restricted                  # Requested from OECD Regpat.
├── codes/
│   ├── 01_network_construction.py  # Network of Places transformation & Graph building
│   ├── 02_novelty_measures.py      # Computation of Uzzi-type Z-scores for atypicality
│   ├── 03_brokerage_metrics.py     # Calculation of Burt's Constraint & Betweenness Centrality
│   ├── 04_panel_regression.R       # TWFE Econometric models & Interaction analysis
│   └── 05_visualization.py         # Generation of maps and correlation plots
├── figures/
│   ├── fig1_zscore_dist.png        # Distribution of IPC pair Z-scores
│   ├── fig2_broker_geography.png   # Map of low-constraint inventor shares
│   └── fig4_interaction.png        # Scatter plot of Brokerage vs. Atypicality
├── documents/
│   └── Virág_Bruinsma_Paper.pdf    # Final research paper
└── README.md
```

### Scripts
The analysis is divided into sequential notebooks corresponding to the methodology sections of the paper.

**file: `01_network_construction.py`**
-   Constructs the co-inventorship networks for each NUTS 2 region using OECD REGPAT data.
-   Applies the **"Network of Places"** transformation (Pizarro, 2007) to mitigate spurious density from large teams.
-   Collapses structurally equivalent inventors (cliques) into single "place" nodes before calculating metrics.

**file: `02_novelty_measures.py`**
-   Calculates the "atypicality" of technology pairs using the method established by Uzzi et al. (2013).
-   Constructs a null model based on observed IPC code frequencies to calculate Z-scores for every pair.
-   Identifies atypical patents (containing pairs with $Z < 0$) and aggregates them to the regional level ($Y_{rt}$).

**file: `03_brokerage_metrics.py`**
-   Calculates **Burt’s Constraint** to identify "Share of Brokers" (inventors in the bottom 25% of the constraint distribution).
-   Calculates **Mean Betweenness Centrality** (log-transformed) to measure the extent to which regional inventors lie on shortest paths.
-   Computes regional network controls: Density, Isolate Share, and Community Structure (Louvain modularity).

**file: `04_econometric_analysis.R`**
-   Implements **Two-Way Fixed Effects (TWFE)** regression models to predict atypical patenting.
-   Includes Region ($\mu_r$) and Time-Window ($\lambda_t$) fixed effects to control for unobserved heterogeneity and common shocks.
-   Tests the interaction between **Brokerage** and **Relatedness Density** to examine if coherent knowledge bases amplify the effect of brokers.

### Data
The project utilizes data derived from the **[OECD REGPAT database](https://www.wipo.int/en/web/economics/research)** (January 2024 edition).
-   Contains patent applications filed with the EPO, restricted to EU and EFTA countries.
-   **Temporal Scope:** 1979–2023, divided into nine non-overlapping 5-year windows.
-   **Spatial Unit:** NUTS 2 regions for analysis; NUTS 3 inventor locations used for granular network construction.

### Figures and Reports

**file: `figures/fig1_zscore_dist.png`**
-   [cite_start]Histogram showing the distribution of IPC pair Z-scores, illustrating the "fat tail" of atypical combinations (negative Z-scores).

**file: `figures/fig2_broker_geography.png`**
-   Choropleth map visualizing the "Geography of Brokerage," specifically the share of low-constraint inventors across European NUTS 2 regions.

**file: `figures/fig4_interaction_scatter.png`**
-   Scatter plot visualizing the relationship between the share of brokers and the share of atypical patents, highlighting regions like DE71 and DE12.

**file: `documents/Virág_Bruinsma_Paper.pdf`**
-   The final research paper detailing the null results for short-term flows and the significant positive findings for long-term novelty capacity.

### Licence
MIT License (MIT): see the [License File](https://github.com/sensiolabs/GotenbergBundle/blob/1.x/LICENSE) for more details.

