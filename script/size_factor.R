#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

### data dir
work.dir <- args[1]
samplesheet <- args[2]
gtf <- args[3]


### test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least one argument is not supplied (input file).n", call.=FALSE)
} 

### library
library(dplyr)
library(tidyr)
library(magrittr)
library(data.table)
library(ggplot2)
library(DESeq2)
library(rtracklayer)
library(GenomicFeatures)

### dir
featureCounts.dir <- paste0(work.dir, "featureCounts/")
data.dir <- paste0(work.dir, "analysis/")

### samplesheet
sample.info <- read.table(samplesheet, header = TRUE)

## Annotation of gene id and gene names
gtf.dat  <- rtracklayer::import(gtf)
gtf.dat  <- as.data.frame(gtf.dat)
gtf.dat %<>% mutate(species = ifelse(seqnames %in% paste0("chr", c(seq(1, 22), "X", "Y", "M")), "human", "yeast"))
gene.dat <- gtf.dat[gtf.dat$type %in% "gene",]
cat("Number of annotated genes in each transcriptome: \n")
table(gene.dat$species)

### count data
raw.data <- list()
count.data <- list()

for (sample in sample.info$sample_name){
  print(sample)
  raw.data[[sample]] <- read.table(file.path(featureCounts.dir, paste0(sample, ".gene.featureCounts.txt")), header = T)
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
count.data <- count.data[,-1]
cat("Number of genes in count data:", nrow(count.data), "\n")

### filter genes that do not express in any samples
count.data <- count.data[apply(count.data, 1, sum) > 0,]
cat("Number of genes in count dats after filtering genes with 0 counts across samples:", nrow(count.data), "\n")

### spike-in ratio 
spikein.count.data <- count.data[rownames(count.data) %in% gene.dat$gene_id[gene.dat$species!="human"], ]
human.count.data <- count.data[rownames(count.data) %in% gene.dat$gene_id[gene.dat$species=="human"], ]
count.ratio <-
  data.frame(sample_id = sample.info$sample_id,
             spikein = colSums(spikein.count.data) / colSums(count.data))

sample.info %<>% left_join(count.ratio)

### DESeq2
#### size factor 
meta.data <- sample.info
rownames(meta.data) <- sample.info$sample_id
dds <- 
  DESeqDataSetFromMatrix(countData = as.matrix(count.data),
                         colData = sample.info,
                         design = ~ 0 + sample)

### normalization to the yeast spike-ins 
dds <- estimateSizeFactors(dds, controlGenes = rownames(count.data) %in% rownames(spikein.count.data))
cat("Normalization factors for each sample: \n")
yeast_size_factor <- data.frame(sample_id = names(sizeFactors(dds)), size_factor=sizeFactors(dds))
yeast_size_factor
write.table(yeast_size_factor, paste0(data.dir, "size_factors_deseq2.txt"), quote=F, sep="\t", row.names =FALSE)
