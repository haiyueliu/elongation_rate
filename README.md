# Gene-specific transcription elongation rates

This repository contains the main script to calculate the absolute RNA Pol II elongation rates for individule genes based on DRB/TTchem-seq2 data.
[DRB/TTchem-seq2]() is an updated version of DRB/TTchem-seq [Gregersen et al., 2020](DOI: 10.1038/s41596-019-0262-3) data. It differs from DRB/TTchem-seq in two main manners:

1. The DRB release times have been reduce from 10, 20, 30 and 40 minutes to 5, 10, 15 and 20 minutes.
2. In stead of calling the peaks of elongation waves, the wave fronts have been called and represent the distance travelled by RNAPII.

We have shown that we increased the number of estimation for individual genes by ~8 fold comparing to DRB/TTchem-seq. The bash and R scripts for data processing and elongation rates calculation are shown below.

## I. Data processing (bash scripts)

### Requirements

### How to start

Download the data_processing.sh script modify the configuration parameters to your own system. The run the script in a server with at least 32G memory.


## II. Elongation rates calculation (R scripts)

### Requirements

Install the required packages:
     rtracklayer
     GenomicFeatures
     tibble
     dplyr
     tidyr
     magrittr
     ggplot2
     cowplot

### How to start 

1.  Download the help_functions.R in the R folder in this repository and source it.
2.  Follow the instructions in elongation_rates.Rmd.
 

## Reference 


