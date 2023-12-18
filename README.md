# Gene-specific transcription elongation rates

This repository contains the main scripts to calculate the absolute RNA Pol II elongation rates for individual genes based on DRB/TTchem-seq2 data.
DRB/TTchem-seq2[1] is an updated version of DRB/TTchem-seq[2], which largely (6-10 folds) increases the number of individual genes whose elongation rates are calculated. The v2 mainly differs from the original version in two aspects:

1. The DRB release times have been reduce from 10, 20, 30 and 40 minutes to 5, 10, 15 and 20 minutes.
2. In stead of calling the peaks of elongation waves, the wave fronts have been called and represent the distance traveled by RNAPII.

To obtain the elongation rates for singles genes from raw DRB/TTchem-seq2 data, we process the raw sequencing reads using a bash script and then calculated the elongation rates using R scrips. 

## I. Data processing (bash scripts)

### Requirements

```
umi_tools
samtools
trimgalore
cutadapt
fastqc
STAR
subread
bedtools
bedGraphToBigWig
R (DESeq2; rtracklayer; GenomicFeatures; dplyr; tidyr; magrittr; data.table; ggplot2)
python (pysam)
```

### How to start

Download the DRB_TTchem_seq_data_processing.sh script in the script folder and modify the configuration parameters. Then run the script in UNIX system with at least 32G memory.

### Main procedures

```
UMI attach (optional)
fastQC & adaptor trimming
STAR alignment
Extract uniquely mapping reads
Deduplication (optional)
Removal of exon-intron-exon mapping reads
Split forward and reverse reads
Count mapped reads
Quantify gene count
Calculate spike-in size factors
Convert bam to bigwig
```


## II. Elongation rates calculation (R scripts)

### Requirements

Install the required R packages:

```
rtracklayer
GenomicFeatures
tibble
dplyr
tidyr
magrittr
ggplot2
cowplot
```

### How to start 

1.  Download the help_functions.R in the R folder and source it.
2.  Follow the instructions in elongation_rates.Rmd to calculate gene elongation rates.


## Reference 

1. 
2. Gregersen LH, Mitter R & Svejstrup JQ (2020) [Using TTchem-seq for profiling nascent transcription and measuring transcript elongation.](https://doi.org/10.1038/s41596-019-0262-3) Nat Protoc 15: 604–627
