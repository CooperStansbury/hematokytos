# Alternative Splicing Events (hg38)

This package provides a comprehensive collection of alternative splicing events for the human genome (hg38). These events have been compiled from various annotation files used by popular alternative splicing quantification tools.

## Description

The `alternativeSplicingEvents.hg38` package contains a data frame of alternative splicing events for humans. The data is sourced from annotation files used by the following tools:

* **MISO**
* **VAST-TOOLS**
* **SUPPA**
* **rMATS**

This resource is useful for researchers working on alternative splicing analysis in human samples.

## Installation

To use this package, you need to have **R** (version "4.5" or higher) and **Bioconductor** installed.

1.  **Install BiocManager:**
    ```R
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    ```

2.  **Install the package:**
    ```R
    BiocManager::install("alternativeSplicingEvents.hg38")
    ```

## Usage

After installation, you can load the data into your R session like any other library.

```R
library(alternativeSplicingEvents.hg38)
```

## Source

[Bioconductor: alternativeSplicingEvents.hg38](https://bioconductor.org/packages/release/data/annotation/html/alternativeSplicingEvents.hg38.html)

## Environment Setup 🛠️

To ensure a clean and reproducible workflow, it's best to create a dedicated Conda environment for this analysis. This environment, named `altsplicing`, will contain the required Bioconductor package and all its dependencies.

Use the following command in your terminal to create and populate the environment using **Mamba**:

```bash
mamba create -n altsplicing -c conda-forge -c bioconda bioconductor-alternativesplicingevents.hg38
mamba activate altsplicing
```

After activating the environment you can import the package within R as
shown above.
