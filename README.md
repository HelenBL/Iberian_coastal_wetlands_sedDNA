# Long-term biodiversity responses across Iberian coastal wetlands

This repository contains the **data and analysis scripts** associated with the research article:

> **Long-term biodiversity responses to anthropogenic change are spatially and taxonomically heterogeneous across Iberian coastal wetlands**

## Overview

Coastal wetlands have undergone centuries of anthropogenic modification, but their long-term biodiversity trajectories remain poorly understood because ecological monitoring generally covers only recent decades.

In this study, we use sedimentary DNA (*seda*DNA) to reconstruct long-term ecological change across Atlantic and Mediterranean coastal wetlands of the Iberian Peninsula. We combine **COI and 18S *seda*DNA metabarcoding** with **<sup>210</sup>Pb and radiocarbon dating, X-ray fluorescence (XRF) geochemistry, and stable-isotope analyses (δ<sup>13</sup>C and δ<sup>15</sup>N)** across eight sediment records.

The analyses examine relationships between eukaryotic biodiversity and anthropogenic enrichment using a general **Anthropogenic Enrichment Index (AEI)**, together with temporal trends, community differentiation, taxonomic composition, and concurrent changes among biological, isotopic, and geochemical proxies.

Overall, the results show that long-term biodiversity responses to anthropogenic enrichment are spatially, temporally, and taxonomically heterogeneous, highlighting the importance of local environmental context and the value of multiproxy palaeoecological approaches for establishing ecological baselines.

## Repository structure

```text
.
├── Data/
│   └── Datasets used in the analyses presented in the manuscript
│
├── Scripts/
│   └── R scripts used for data processing, statistical analyses, and figure generation
│
├── README.md
└── LICENSE
```

## Data

The `Data/` directory contains the datasets used in the analyses presented in the manuscript and its supplementary material.

The available datasets include:

- **sedaDNA metabarcoding data (COI):** processed sedimentary DNA metabarcoding data for the COI genetic marker.
- **sedaDNA metabarcoding data (18S):** processed sedimentary DNA metabarcoding data for the 18S genetic marker.
- **Sediment chronology data:** <sup>210</sup>Pb and radiocarbon (<sup>14</sup>C) dating data used to establish the age–depth models of the sediment records.
- **Anthropogenic Enrichment Index (AEI):** data used to quantify relative multielement enrichment above site- and element-specific historical baselines.
- **X-ray fluorescence (XRF) data:** geochemical data used to characterise temporal changes in elemental composition and anthropogenic enrichment.
- **Stable-isotope data:** δ<sup>13</sup>C and δ<sup>15</sup>N measurements used as complementary environmental proxies in the multiproxy reconstructions.

Together, these datasets provide the biological, chronological, geochemical, and isotopic information required to reproduce the analyses, tables, and figures presented in the article.

## Scripts

The `Scripts/` directory contains the R scripts used for data processing, statistical analyses, and figure generation.

The scripts are organised as follows:

# `01_Age_depth_models_210Pb_14C.R`

Builds the sediment age–depth models by combining ^210Pb and AMS ^14C dating information using `rplum`.

# `02_XRF_Anthropogenic_Enrichment_Index.R`

Processes the XRF geochemical data and performs the analyses used to characterise anthropogenic elemental enrichment and calculate the Anthropogenic Enrichment Index (AEI).

# `03_iNEXT_sampling_completeness.R`

Performs iNEXT-based analyses to evaluate sequencing/sampling completeness and diversity patterns for the COI and 18S metabarcoding datasets.

# `04_Richness_AEI_temporal_analysis.R`

Calculates eukaryotic ASV richness for the COI and 18S markers and evaluates relationships between richness and the Anthropogenic Enrichment Index (AEI), including temporal trends and rolling correlations.

This script generates the analyses and outputs associated with **Figure 1**.

# `05_Taxonomic_richness_through_time.R`

Examines temporal changes in taxonomic richness and composition for the COI and 18S markers, including total eukaryotic richness and the dominant taxonomic groups.

This script generates the analyses and outputs associated with **Figure 2**.

# `06_nMDS_PERMANOVA_community_composition.R`

Analyses differences in community composition using non-metric multidimensional scaling (nMDS) and PERMANOVA for COI and 18S datasets.

Analyses include all retained ASVs and taxonomic subsets such as Eukaryota and Metazoa.

This script generates the analyses and outputs associated with **Figure 3**.

# `07_Multiproxy_temporal_reconstruction.R`

Integrates COI and 18S ASV richness, stable-isotope records, and the Anthropogenic Enrichment Index (AEI) to generate multiproxy temporal reconstructions across sediment records.

This script generates the analyses and outputs associated with **Figure 4**.


## Reproducibility

The scripts are numbered according to their logical position in the analytical workflow.

Analyses can be reproduced using the datasets available in the `Data/` directory and the R scripts provided in `Scripts/`.

Software versions, R package dependencies, and additional instructions required to reproduce individual analyses can be added to this section.
