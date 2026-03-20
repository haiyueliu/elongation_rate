#!/usr/bin/env Rscript
#############################
### Helper Functions
#############################
### Load libraries
suppressMessages(suppressWarnings(library(rtracklayer)))
suppressMessages(suppressWarnings(library(GenomicFeatures)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(magrittr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))

### Define helper operator
`%notin%` <- Negate(`%in%`)

##########################################
### 1. Filter overlapping gene annotations
##########################################
filter_overlapping_genes <- function(gtf          = gtf.file,        ### full path to annotation gtf (Gencode / Ensembl)
                                     lower.length = min.genewidth,
                                     upper.length = max.genewidth,
                                     count.data   = NULL) {

  ### Read gtf as GRanges
  gtf <- rtracklayer::import(gtf)
  gene.gr <- gtf[gtf$type == "gene"]
  gene.gr$gene_width <- width(gene.gr)

  ### Filter genes that overlap other genes on the same strand.
  ### If count.data is provided, restrict to expressed genes first.
  if (!is.null(count.data)) {
    count.data <- count.data[rowSums(count.data) >= 2, ]
    gene.gr.expressed <- gene.gr[gene.gr$gene_id %in% rownames(count.data)]
    gene.gr <- gene.gr.expressed[countOverlaps(gene.gr.expressed, gene.gr.expressed, ignore.strand = FALSE) == 1, ]
  } else {
    gene.gr <- gene.gr[countOverlaps(gene.gr, gene.gr, ignore.strand = FALSE) == 1, ]
  }

  ### Keep genes within specified length range
  gene.gr <- gene.gr[width(gene.gr) >= lower.length & width(gene.gr) <= upper.length]

  return(gene.gr)
}

#######################################################
### 2. Filter lowly expressed genes
###    Compares mean coverage per gene across samples
#######################################################
# FIX: added sample.info as an explicit parameter instead of relying on global scope
filter_low_expr_genes <- function(cov.list        = cov.list,
                                  mean.cov.cutoff = 0.01,
                                  control.samples = control.sample.names,
                                  sample.info     = sample.info) {

  ### Calculate mean coverage across all positions for each gene
  mean.cov.list <- lapply(cov.list, function(x) { unlist(lapply(x, function(y) { mean(y) })) })
  mean.cov.df   <- as.data.frame(mean.cov.list)

  mean.cov.df %<>%
    rownames_to_column("gene_id") %>%
    pivot_longer(cols = -gene_id, names_to = "sample", values_to = "mean.cov") %>%
    mutate(sample = factor(sample, levels = sample.info$sample_name))  # now uses explicit parameter

  ### Distribution plot
  p_mean_cov <-
    mean.cov.df %>%
    ggplot(aes(mean.cov)) +
      geom_histogram(bins = 50, color = "grey") +
      geom_vline(xintercept = mean.cov.cutoff, linetype = "dashed") +
      scale_x_log10() +
      scale_y_log10() +
      facet_wrap(~sample, ncol = 2) +
      theme_bw()

  ### Keep genes that pass the coverage cutoff in all non-control samples
  n_noncontrol <- length(setdiff(names(cov.list), control.samples))
  genes.keep <-
    mean.cov.df %>%
    filter(sample %notin% control.samples) %>%
    filter(mean.cov >= mean.cov.cutoff) %>%
    group_by(gene_id) %>%
    summarise(n_sample = n()) %>%
    ungroup() %>%
    filter(n_sample == n_noncontrol)   # FIX: use computed count, not length(cov.list) - length(control.samples)

  genes.keep <- genes.keep$gene_id

  return(list(
    mean.cov.df   = mean.cov.df,
    mean.cov.plot = p_mean_cov,
    genes.keep    = genes.keep
  ))
}

##################################
### 3. Wave front calling
###    Input: smoothed spline coverage vector for one gene (TSS -> TES)
##################################
# FIX: removed unused parameter min.cosine
wave_front_calling <- function(y,
                               offset = 1000) {   ### wave fronts must be >= 1 kb downstream of TSS

  ### Cumulative sum of coverage
  y <- cumsum(y)

  ### Normalize x and y to the unit square [0, 1]
  y_len <- length(y)
  x <- seq(y_len) / y_len
  y <- y / max(y)

  ### Euclidean distance between first and last point on the cumulative curve
  dmax <- sqrt((x[1] - x[y_len])^2 + (y[1] - y[y_len])^2)

  ### Distances from each point to the first and last points
  d1 <- sqrt((x - x[1])^2    + (y - y[1])^2)
  d2 <- sqrt((x[y_len] - x)^2 + (y[y_len] - y)^2)

  ### Cosine similarity via the law of cosines
  # FIX: guard against division by zero at endpoints (d1=0 or d2=0 -> NaN)
  cosine_values <- abs((d1^2 + d2^2 - dmax^2) / (2 * pmax(d1 * d2, 1e-10)))

  ### Index of minimum cosine (= sharpest bend = wave front),
  ### constrained to be at least `offset` positions downstream of TSS
  whichmin <- offset + which.min(cosine_values[seq(offset + 1, y_len)])

  ### Return 0 if no valid wave front found
  whichmin <- ifelse(is.na(whichmin),       0, whichmin)
  whichmin <- ifelse(is.logical(whichmin),  0, whichmin)

  return(whichmin)
}
