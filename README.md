# Do Regions with More Brokers Produce More Atypical Technology Combinations?
## Investigating the Role of Inventor Brokerage and Small-World Networks in Regional Innovation

### Course: Economic Complexity – Social Data Science (MSc)

### Team Members:
-   Patrik Bruinsma
-   Ákos Virág

### Project Overview
This repository contains the research pipeline for a project analyzing the relationship between network structure and regional innovation novelty. Innovation often emerges from the recombination of existing technologies, yet only a small fraction of inventions relies on atypical combinations of knowledge domains.

The primary objective of this research is to bridge firm-level brokerage mechanisms with regional-level outcomes. We investigate whether regions with a higher density of "inventor brokers" (individuals bridging structural holes) in their co-inventor networks generate a higher share of atypical technology combinations in subsequent years. Furthermore, we examine if this effect is moderated by the region's "small-world" network structure.

### Repository Structure
```text
.
├── data/
│   ├── patstat_sample_eu27.csv        # (Restricted) Subset of EPO patent data for EU27+4 regions 
│   └── regional_ipc_matrix.csv        # Processed Region-Year-Technology matrix
├── codes/
│   ├── 01_network_construction.ipynb  # Construction of co-inventor networks
│   ├── 02_brokerage_metrics.ipynb     # Calculation of Burt's constraint and broker density
│   ├── 03_novelty_measures.ipynb      # Computation of Uzzi-type z-scores for atypicality
│   ├── 04_panel_regression.ipynb      # Main econometric analysis
│   └── 05_main_pipeline.py            # Execution script for the full workflow 
├── figures/
│   ├── small_world_interaction.png    # Interaction plots of brokerage and network structure
│   └── atypicality_trends.png         # Visualization of atypical patent shares over time
├── documents/
```

### Scripts
The analysis is divided into sequential notebooks corresponding to the methodology sections of the paper.

**file: `01_network_build.ipynb`**
-   Constructs the cumulative co-inventor network from the OECD REGPAT database.
-   Filters for the "active" inventor population (Degree $\geq 2$) to exclude isolates.
-   Calculates Burt’s Constraint (Eq. 1) to identify Weak ($C < \mu$) and Strong ($C \leq \mu - \sigma$) brokers.

**file: `02_inequality.ipynb`**
-   Calculates the Herfindahl-Hirschman Index (HHI) to measure the concentration of brokerage power across 31 European countries and IPC4 technological fields.
-   Generates Lorenz curves to visualize distributional inequality, finding high geographic concentration (HHI $\approx 0.79$ for strong brokers).

**file: `03_robustness.ipynb`**
-   Simulates reverse percolation-style targeted node removal on the giant connected component (GCC)[cite: 63].
-   Compares the impact of removing top-degree brokers ($k=100, 500, \dots, 10000$) versus random nodes on the average shortest path length to assess network fragility[cite: 66, 67, 68].
-   Targeted removal increases path length by $\approx 40\%$, whereas random removal shows minimal impact[cite: 70, 72].

**file: `04_community.ipynb`**
-   Applies the Louvain algorithm to detect community structures, finding a high modularity of $0.994$.
-   Determines the core-periphery placement of brokers using k-core decomposition.
-   Definitions used: Periphery ($k_{core} \leq 2$), Core ($k_{core} \geq 5$), Intermediate (otherwise).

**file: `05_broker_typology.ipynb`**
-   Implements feature-based role discovery using K-means clustering ($k=4$) on strong brokers.
-   Identifies four broker profiles:
    -   **Type 1:** High clustering (embedded in cohesive teams, mean local clustering $\approx 0.95$).
    -   **Type 2:** High degree (super-star hubs, mean degree $\approx 125.5$).
    -   **Type 3:** Low constraint (bridging connectors, constraint $\approx 0.16$).
    -   **Type 4:** Other (peripheral/weakly embedded).

### Data
The project utilizes data derived from the **OECD REGPAT database** (January 2024 edition).

**file: `regpat_nodes.csv` & `regpat_edges.csv`**
-   Represents the co-inventor network where nodes are inventors and edges are weighted by shared patents.
-   The analysis is restricted to inventors with at least two co-inventor links (Degree $\ge 2$) to exclude incidental one-time inventors.

### Figures and Reports

**file: `figures/backbone_map.png`**
-   Visualizes the maximum spanning tree of NUTS2-level collaboration ties, highlighting high-activity regions like Germany and France (Figure 1).

**file: `figures/broker_geography.png`**
-   Maps the clustering of countries based on their broker role composition, revealing distinct groups (e.g., Western/Northern economies vs. Eastern Europe) (Figure 2).

**file: `documents/network_groupwork.pdf`**
-   The final research paper detailing the finding that European brokers form a dense "Elite Club" and that brokerage is technologically specialized rather than a generalist phenomenon.

```tex
.
└── research_plan.pdf              # Full research proposal and literature review
├── requirements.txt                   # Python dependencies (networkx, pandas, statsmodels)
└── README.md                          # Project description and usage
```

### Licence
MIT License (MIT): see the [License File](https://github.com/sensiolabs/GotenbergBundle/blob/1.x/LICENSE) for more details.

