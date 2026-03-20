#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

### Arguments
if (length(args) < 3) {
  stop("Not all arguments are supplied. Usage: size_factor.R <work.dir> <samplesheet> <gtf>\n", call. = FALSE)
}
work.dir   <- args[1]
samplesheet <- args[2]
gtf        <- args[3]

### Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(magrittr)
  library(data.table)
  library(ggplot2)
  library(DESeq2)
  library(rtracklayer)
  library(GenomicFeatures)
})

### Directories
featureCounts.dir <- file.path(work.dir, "featureCounts")
data.dir          <- file.path(work.dir, "analysis")

### Samplesheet
# FIX: typo 'atringAsFactors' -> 'stringsAsFactors'
sample.info <- read.table(samplesheet, header = TRUE, stringsAsFactors = FALSE)

### Annotation
gtf.dat <- rtracklayer::import(gtf)
gtf.dat <- as.data.frame(gtf.dat)
gtf.dat %<>% mutate(species = ifelse(seqnames %in% paste0("chr", c(seq(1, 22), "X", "Y", "M")), "human", "yeast"))
gene.dat <- gtf.dat[gtf.dat$type %in% "gene", ]
cat("Number of annotated genes by species:\n")
print(table(gene.dat$species))

### Read count data
raw.data   <- list()
count.data <- list()
for (sample in sample.info$sample_name) {
  print(sample)
  # FIX: removed trailing comma inside file.path(); fixed typo 'stringAsFactors' -> 'stringsAsFactors'
  raw.data[[sample]] <- read.table(
    file.path(featureCounts.dir, paste0(sample, ".gene.featureCounts.txt")),
    header = TRUE, stringsAsFactors = FALSE
  )
  colnames(raw.data[[sample]]) <- c("gene_id", "chr", "start", "end", "strand", "length", "count")
  count.data[[sample]] <-
    raw.data[[sample]] %>%
    dplyr::select(gene_id, count) %>%
    mutate(sample = sample)
}

count.data <- rbindlist(count.data)
count.data %<>%
  pivot_wider(names_from = sample, values_from = count) %>%
  as.data.frame()
rownames(count.data) <- count.data$gene_id
count.data <- count.data[, -1]

cat("Number of genes in count data:", nrow(count.data), "\n")

### Filter genes with 0 counts across all samples
count.data <- count.data[apply(count.data, 1, sum) > 0, ]
cat("Number of genes after filtering 0-count genes:", nrow(count.data), "\n")

### Spike-in vs human split
spikein.count.data <- count.data[rownames(count.data) %in% gene.dat$gene_id[gene.dat$species != "human"], ]
human.count.data   <- count.data[rownames(count.data) %in% gene.dat$gene_id[gene.dat$species == "human"], ]

### Spike-in ratio (informational)
count.ratio <-
  data.frame(
    # FIX: use sample_name (matches count.data colnames) so the left_join key aligns
    sample_name = sample.info$sample_name,
    spikein     = colSums(spikein.count.data) / colSums(count.data)
  )
sample.info %<>% left_join(count.ratio, by = "sample_name")

### DESeq2 size factors
# FIX: rownames of colData must match colnames of countData (both are sample_name)
sample.info$sample <- as.factor(sample.info$sample)   # FIX: design variable must be a factor
rownames(sample.info) <- sample.info$sample_name

dds <-
  DESeqDataSetFromMatrix(
    countData = as.matrix(count.data),
    colData   = sample.info,
    design    = ~ 0 + sample
  )

### Normalize to yeast spike-ins
dds <- estimateSizeFactors(dds, controlGenes = rownames(count.data) %in% rownames(spikein.count.data))

cat("Normalization factors for each sample:\n")
# FIX: use sample_name as column 1 (bash step 11 awk matches on sample_name, not sample_id)
yeast_size_factor <- data.frame(
  sample_name = names(sizeFactors(dds)),
  size_factor = sizeFactors(dds)
)
print(yeast_size_factor)

# FIX: use file.path() to avoid missing separator (paste0 would give .../analysissize_factors...)
write.table(
  yeast_size_factor,
  file.path(data.dir, "size_factors_deseq2.txt"),
  quote     = FALSE,
  sep       = "\t",
  row.names = FALSE
)
