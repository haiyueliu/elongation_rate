# Gene-specific transcription elongation rates

This repository contains the main scripts to calculate the absolute RNA Pol II elongation rates for individule genes based on DRB/TTchem-seq2 data.
DRB/TTchem-seq2 [1] is an updated version of DRB/TTchem-seq (Gregersen et al., 2020)[2] data. It mainly differs from DRB/TTchem-seq in two folds:

1. The DRB release times have been reduced from 10, 20, 30 and 40 minutes to 5, 10, 15 and 20 minutes.
2. In stead of calling the peaks of elongation waves, the wave fronts have been called and represent the distance travelled by RNAPII.

We have shown that the number of estimation for individual genes has increased by ~8 fold (i.e. from hundreds to thousands) comparing to DRB/TTchem-seq. Follow the bash and R scripts for data processing and elongation rates calculation as shown below.

## I. Data processing (bash scripts)

### Requirements


### How to start

Download the data_processing.sh script and modify the configuration parameters to your own system. 
Then run the bash script in a server with at least 32G memory.


## II. Elongation rates calculation (R scripts)

### Requirements

Install the required R packages:
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
2.  Follow the instructions in elongation_rates.Rmd to calculate gene elongation rates.
 

## Reference 
1.
2. Gregersen LH, Mitter R & Svejstrup JQ (2020) [Using TTchem-seq for profiling nascent transcription and measuring transcript elongation.](https://doi.org/10.1038/s41596-019-0262-3) Nat Protoc 15: 604–627
2. 

