#!/usr/bin/env Rscript
#############################
### Helper Functions
#############################
### load libraries
suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(GenomicFeatures)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
### define function 
`%notin%` <- Negate(`%in%`)

##########################################
### 1. function to filter gene annotation 
##########################################
filter_overlapping_genes <- function(gtf = gtf.file,                    ### full path to annotation file in gtf format (Gencode / Ensembl)
                                     lower.length = min.genewidth,
                                     upper.length = max.genewidth,
                                     count.data = NULL) {
  
  ### read gtf file as gene range file
  gtf <- rtracklayer::import(gtf.file)
  ### only keep gene tracks
  gene.gr <- gtf[gtf$type=="gene"]
  gene.gr$gene_width <- width(gene.gr)
  
  ### filter genes that overlap with other genes on the same strand. If count_data if provided, only look at genes that are expressed in the data
  if (!is.null(count.data)) {
    count.data <- count.data[rowSums(count.data) >= 2, ]
    gene.gr.expressed <- gene.gr[gene.gr$gene_id %in% rownames(count.data)]
    gene.gr <- gene.gr.expressed[countOverlaps(gene.gr.expressed, gene.gr.expressed, ignore.strand=FALSE) == 1, ]
  } else{
    gene.gr <- gene.gr[countOverlaps(gene.gr, gene.gr, ignore.strand=FALSE) == 1, ]
  }
  
  ### keep genes with certain lengths 
  gene.gr <- gene.gr[width(gene.gr) >= lower.length & width(gene.gr) <= upper.length]
  
  ### return gene range object
  return(gene.gr)
}

#######################################################
### 2. filter lowly expressed gene function 
### compare the mean coverage per gene across samples 
#######################################################
filter_low_expr_genes <- function(cov.list = cov.list,
                                  mean.cov.cutoff = 0.01,
                                  control.samples = control.sample.names){
  ### calculate the mean coverage across all position for each gene
  mean.cov.list <- lapply(cov.list, function(x){ unlist(lapply(x, function(y) {mean(y)})) })
  mean.cov.df <- as.data.frame(mean.cov.list)
  mean.cov.df %<>%
    rownames_to_column("gene_id") %>%
    pivot_longer(cols=-gene_id, names_to = "sample", values_to = "mean.cov") %>%
    mutate(sample = factor(sample, levels = sample.info$sample_name))
  ### distributions
  p_mean_cov <-
    mean.cov.df %>%
    ggplot(aes(mean.cov)) %>%
    + geom_histogram(bins=50, color="grey") %>%
    + geom_vline(xintercept = mean.cov.cutoff, linetype="dashed") %>%
    + scale_x_log10() %>%
    + scale_y_log10() %>%
    + facet_wrap(~sample, ncol=2) %>%
    + theme_bw()
  ### keep highly expressed genes
  genes.keep <- 
    mean.cov.df %>%
    filter(sample %notin% control.samples) %>%
    filter(mean.cov >= mean.cov.cutoff) %>%
    group_by(gene_id) %>%
    summarise(n_sample = n()) %>%
    ungroup %>%
    filter(n_sample == length(cov.list) - length(control.samples))
  genes.keep <- genes.keep$gene_id
  ### return kept gene ids
  return(list(mean.cov.df = mean.cov.df,
              mean.cov.plot = p_mean_cov, 
              genes.keep = genes.keep))
}
##################################
### 3. wave front calling
###################################
### input y is a vector: smoothed spline of normalized coverage for each gene between TSS and TES
wave_front_calling <- function(y,               
                               min.cosine = 1, 
                               offset = 1000){  ### add one constraint (wave fronts need to be at least 1kb downstream of the TSS)
  
  
  ### Calculate the cumulative sums of the input vector
  y <- cumsum(y)
  
  ### Normalize x and y values to the unit square (0 to 1).
  y_len<- length(y)
  x <- seq(y_len) / y_len
  y <- y / max(y)
  
  ### Calculate the Euclidean distance between the first point and the last point on the cumulative curve
  dmax <- sqrt((x[1] - x[y_len])^2 + (y[1] - y[y_len])^2)
  
  ### Calculate distances to the first and last points
  d1 <- sqrt((x - x[1])^2 + (y - y[1])^2)
  d2 <- sqrt((x[y_len] - x)^2 + (y[y_len] - y)^2)
  
  ### Calculate cosine similarity (using vectorization to avoid looping through every points in the cumulative curve)
  cosine_values <- abs((d1^2 + d2^2 - dmax^2) / (2 * d1 * d2))
  
  ### Find the index of the minimum cosine value
  ### add one constraint (wave fronts need to be at least 1kb downstream of the TSS)
  whichmin <- offset + which.min(cosine_values[seq(offset+1, y_len)])
  ### whichmin <- ifelse(is.nae(whichmin), 0, whichmin)
  return(whichmin)
}
