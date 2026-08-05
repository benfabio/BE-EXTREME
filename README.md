<img width="3780" height="947" alt="EUROPE-Germany-Berlin-1850-2025-BK" src="https://github.com/user-attachments/assets/a2715b24-6c18-4847-a972-b70cd9078b36" />

# BE-EXTREME

Here are the R scripts developed to build the **BE-EXTREME** dataset and reproduce the analyses and figures presented in:

> Benedetti et al. *BE-EXTREME : Climate and extreme event reconstructions (1950-2024) for the grasslands and forests of the Biodiversity Exploratories*. in preparation for *Earth System Science Data* (ESSD).

---

## Overview

This repository contains the main R scripts for building, processing, modelling, and visualizing the data underlying BE-EXTREME. It accompanies (but does not duplicate) the dataset itself, which is archived here **[Add Zenodo link here]**.

This is a curated public subset of a larger private development repository. Scripts here cover data processing, index computation, statistical modelling, and figure generation relevant to the published paper; exploratory and unpublished analyses are not included.

---

Scripts are **NOT** numbered to reflect the intended execution order yet. Pleasr find here a brief description of each script's purpose:

| Script | Purpose |
|---|---|
| `01_data_import.R` | Loads and merges raw climate and ecological input data |
| `02_quality_control.R` | Flags and filters missing/erroneous records |
| `03_spei_computation.R` | Computes SPEI drought indices from climate variables |
| `...` | `...` |
| `19_figures_main.R` | Generates main text figures |

---

## Requirements

- R (≥ 4.5.X recommended)
- Key packages include: `tidyverse, reshape2, lubridate, lme4, mgcv, SPEI, etc.`

Install dependencies with:

```r
install.packages(c("tidyverse","data.table","reshape2","lubridate","zoo","parallel","qmap","lme4","mgcv","terra","SPEI","viridis","scales"))
```

---

## Usage

1. Clone this repo.
2. Update file paths at the top of each script to point to your local copy of the input data (see **[Add Zenodo link here]**).
3. Run scripts to reproduce the processing pipeline, statistical models, and figures from the paper.

Scripts are documented with inline comments; each corresponds to a specific step, table, or figure referenced in the manuscript (cross-references noted where applicable).

---

## Data availability

The BE-EXTREME dataset is archived at **[repository name, e.g. Zenodo/PANGAEA]**: [DOI link]. This code repository expects the archived dataset as input and does not host the raw or processed data files directly.

---

## Citation

If you use this code, please cite the associated preprint:

> Benedetti . et al. (in prep.). *BE-EXTREME : Climate and extreme event reconstructions (1950-2024) for the grasslands and forests of the Biodiversity Exploratories*. [DOI link]

A `CITATION.cff` file is included for automated citation tools (e.g. GitHub's "Cite this repository" feature).

---

## License

This code is released under the **[MIT / GPL-3.0 / CC-BY-4.0 — choose one]** license. See [`LICENSE`](LICENSE) for details.

---

## Contact

For questions, please contact **Fabio Benedetti** (fabio.benedetti@unibe.ch), or open an issue in this repository.

## Climate stripes
The climate stripes shown above were issued from: https://showyourstripes.info/ 

> Graphics and lead scientist: Ed Hawkins, National Centre for Atmospheric Science, UoR.
> Data: Berkeley Earth, NOAA, UK Met Office, MeteoSwiss, DWD, SMHI, UoR & ZAMG

