#!/usr/bin/env Rscript
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



#################################################################
### Step 1. Filter gene annotation 
### Filter overlapping genes and genes above certain gene length ranges
#################################################################
### call the filter overlapping gene dunction to get the gene ranges with non-overlapping genes
gtf.file <- "~/Documents/GitHub/DRB_TTchem_seq2/annotation/gencode.v43.primary_assembly.basic.annotation.gtf"
min.genewidth <- 30000
max.genewidth <- 300000

gene.gr <- filter_overlapping_genes(gtf = gtf.file,                     ### full path to annotation file in gtf format (Gencode / Ensembl)
                                    lower.length = min.genewidth,
                                    upper.length = max.genewidth,
                                    count_data = NULL)
gene.gr$strand <- strand(gene.gr)
cat("The number of genes left after filtering ovelapping genes and limiting gene lengths:", length(gene.gr))

gene.gr.fwd <- gene.gr[strand(gene.gr) == "+"]
gene.gr.rev <- gene.gr[strand(gene.gr) == "-"]

#######################################################################################
### Step 2. Count single base pair resolution reads depth for every filtered gene
### Count the raw read depth for each gene from its TSS to its TES from bam files ( bamCoverage function from bamsignals package)
### Normalize the read depth to the size factors of each sample 
#######################################################################################
sample.info <- read.table("~/Documents/GitHub/DRB_TTchem_seq2/R/sample.info.txt", header = TRUE)
bam.dir <- "~/Documents/GitHub/DRB_TTchem_seq2/test_data/"
data.dir <- "~/Documents/GitHub/DRB_TTchem_seq2/output/"
control.sample.names <- c("DRB_noR_rep1", "DRB_noR_rep2")
control.sample <- "DRB_noR"

cov.list <- list()
for (sample in sample.info$sample_name){
  print(sample)
  ### retrieve the size factor for each sample from sample.info table
  size.factor = sample.info$size_factor[sample.info$sample_name==sample]
  bam.fwd.file = paste0(bam.dir, sample, ".test.fwd.bam")
  bam.rev.file = paste0(bam.dir, sample, ".test.rev.bam")
  ### count read depth for genes on the forward and reverse strands
  cov.fwd <- bamsignals::bamCoverage(bampath = bam.fwd.file,
                                     gr = gene.gr.fwd,
                                     filteredFlag = 1024,
                                     paired.end = "ignore",
                                     mapqual = 20,
                                     verbose = FALSE)
  cov.rev <- bamsignals::bamCoverage(bampath = bam.rev.file,
                                     gr = gene.gr.rev,
                                     filteredFlag = 1024,
                                     paired.end = "ignore",
                                     mapqual = 20,
                                     verbose = FALSE)
  cov.fwd <- cov.fwd@signals
  cov.rev <- cov.rev@signals
  ### add gene ids to list names
  names(cov.fwd) <- gene.gr.fwd$gene_id
  names(cov.rev) <- gene.gr.rev$gene_id
  ### combine genes on the forward and reverse strands
  cov.list[[sample]] <- c(cov.fwd, cov.rev)
  
  ### remove genes that do not express in these samples
  cov.list[[sample]] <- cov.list[[sample]][lapply(cov.list[[sample]], sum) > 0]
  
  ### normalize read depths to pre-calculated sample size factors (i.e. spike-in size factors or library sizes)
  cov.list[[sample]] <- lapply(cov.list[[sample]], function(x){round(x/size.factor, 3)})
}

######################################################
### Step 3. Remove lowly expressed genes
######################################################
### In our data, we filter genes with mean.cov < 0.01 in any of the non-control samples
### It is optimal that you compare the mean coverage distributions of genes across samples and define the cutoff accordingly.
### Keep genes with certain mean coverage (>=0.01)
genes.keep <- filter_low_expr_genes(cov.list = cov.list, mean.cov.cutoff = 0.01)$genes.keep
for (sample in names(cov.list)){
  cov.list[[sample]] <- cov.list[[sample]][genes.keep]
}
###############################################################################
### Step 4. Normalize coverage values to the control samples
### Subtract the read depth at every base position to the control samples
###############################################################################
#### Mean coverage across replicates
mean.cov.list <- list()
for(sample in unique(sample.info$sample)){
  print(sample)
  ### mean coverage for each position in each gene across replicates
  list1 <- cov.list[[ sample.info$sample_name[sample.info$sample == sample & sample.info$rep == "rep1"] ]]
  list2 <- cov.list[[ sample.info$sample_name[sample.info$sample == sample & sample.info$rep == "rep2"] ]]
  for(m in names(list1)){
    mean.cov.list[[sample]][[m]] <- rowMeans(cbind(list1[[m]], list2[[m]]))
  }
}
norm.cov.list <- list()
### Subtract the mean coverage of each sample by the mean coverage of control sample.
for(sample in names(mean.cov.list)){
  if(sample %notin% control.sample){
    print(sample)
    list.drb <- mean.cov.list[[sample]]
    list.ctrl <- mean.cov.list[[control.sample]]
    for(m in names(list1)){
      norm.cov.list[[sample]][[m]] <- pmax(0, list.drb[[m]] - list.ctrl[[m]])
    }
  }
}
### save data
# saveRDS(norm.cov.list, "control.normalized.cov.list.Rds")
### releases memory
rm(mean.cov.list)
rm(list.drb)
rm(list.ctrl)
gc()

######################################################
### Step 4. Smooth the coverage profiles using splines
#######################################################
gene.spar=0.9     
norm.spline.list <- list()
for (sample in names(norm.cov.list)){
  ### spline for each sample (mean across replicates) after control normalization
  print(sample)
  norm.spline.list[[sample]] <- sapply(norm.cov.list[[sample]], function(x){pmax(0, round(smooth.spline(1:length(x), x, spar=gene.spar)$y, 3))})
}

# ### save data
# saveRDS(norm.cov.list, paste0(data.dir, "per.sample.spike.in.control.normalized.coverage.list.unspliced.Rds"))
# saveRDS(norm.spline.list, paste0(data.dir, "per.sample.spike.in.control.normalized.spline.list.unspliced.Rds"))
# rm(norm.cov.list)
# rm(norm.spline.list)
# gc()
############################################################
### Step 5. Identify wave fronts of the smoothing splines
############################################################
### find the wave front of the elongating RNAPII products relative to the transcript's TSS
wf <- sapply(norm.spline.list, function(x) {sapply(x, function(i) {wave_front_calling(i)})}) %>% as.data.frame
wf[is.na(wf)] <- 0
# write.table(wf, paste0(data.dir, "wf_genes_all.txt"), quote=FALSE, sep="\t", row.names = T, col.names = T)
# wf <- read.table(paste0(data.dir, "wf_genes_all.txt"))

wf_long <-
  wf %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(cols = -gene_id, names_to = "sample", values_to = "wf") %>%
  mutate(sample = factor(sample, levels = unique(sample.info$sample)))


######################################################
### Step 6. Calculate elongation rate for every gene
### filter genes based on wave fronts positions
### linear regression on wave fronts and elongation times
######################################################
# 1. The first wave front should be 0.5kb away from TSS
# 2. At lease three wave fronts should be consecutively increasing.



### order of the wave front position for each gene
wave.orders <- apply(wf, 1, function(x) {paste0(order(x), collapse = "")})

rate_df <-
  wf %>%
  rownames_to_column("gene_id") %>%
  left_join(gene.gr %>%
              as.data.frame %>%
              mutate(gene_length = width) %>%
              dplyr::select(gene_id, gene_name, gene_type, gene_length)) %>%
  pivot_longer(cols = -c(gene_id, gene_type, gene_name, gene_length), names_to = "DRB", values_to = "wf") %>%
  mutate(wf = as.numeric(wf/1000)) %>%
  mutate(time = as.numeric(sub("minR", "", sub("DRB_", "", DRB)))) %>%
  left_join(as.data.frame(wave.orders) %>%
              rownames_to_column("gene_id"))

### recheck the orders of the position of wave fronts
# rate_df %<>%
#   group_by(gene_id) %>%
#   mutate(order = ifelse((wf[DRB=="DRB_5minR"] < wf[DRB=="DRB_10minR"]) &
#                            (wf[DRB=="DRB_10minR"] < wf[DRB=="DRB_15minR"]) &
#                            (wf[DRB=="DRB_15minR"] < wf[DRB=="DRB_20minR"]), "1234",
#                          ifelse(t(wf[DRB=="DRB_5minR"] < wf[DRB=="DRB_10minR"]) &
#                                  (wf[DRB=="DRB_10minR"] < wf[DRB=="DRB_15minR"]), "123",
#                                 ifelse((wf[DRB=="DRB_10minR"] < wf[DRB=="DRB_15minR"]) &
#                                         (wf[DRB=="DRB_15minR"] < wf[DRB=="DRB_20minR"]) &
#                                         (wf[DRB=="DRB_20minR"] != round(gene_length/1000,1)), "234", "other"))))
### recheck the orders of the position of wave fronts
rate_df %<>%
  group_by(gene_id) %>%
  mutate(order = ifelse((wf[DRB=="DRB_5minR"] < wf[DRB=="DRB_10minR"]) &
                          (wf[DRB=="DRB_10minR"] < wf[DRB=="DRB_15minR"]) &
                          (wf[DRB=="DRB_15minR"] < wf[DRB=="DRB_20minR"]) &
                          (wf[DRB=="DRB_20minR"] != round(gene_length/1000, 1)), "1234",
                        ifelse((wf[DRB=="DRB_5minR"] < wf[DRB=="DRB_10minR"]) &
                                 (wf[DRB=="DRB_10minR"] < wf[DRB=="DRB_15minR"]), "123",
                               ifelse((wf[DRB=="DRB_5minR"] < wf[DRB=="DRB_10minR"]) &
                                        (wf[DRB=="DRB_10minR"] < wf[DRB=="DRB_20minR"]), "124",
                                      ifelse((wf[DRB=="DRB_5minR"] < wf[DRB=="DRB_15minR"]) &
                                               (wf[DRB=="DRB_15minR"] < wf[DRB=="DRB_20minR"]), "134",
                                             ifelse((wf[DRB=="DRB_10minR"] < wf[DRB=="DRB_15minR"]) &
                                                      (wf[DRB=="DRB_15minR"] < wf[DRB=="DRB_20minR"]), "234", "other"))))))

rate_df %<>%
  group_by(gene_id) %>%
  mutate(elongation.rate = ifelse(order == "1234", lm(wf~time)$coef[2],
                                  ifelse(order == "123", lm(wf[DRB!="DRB_20minR"] ~ time[DRB!="DRB_20minR"])$coef[2],
                                         ifelse(order == "124", lm(wf[DRB!="DRB_15minR"] ~ time[DRB!="DRB_15minR"])$coef[2],
                                                ifelse(order == "134", lm(wf[DRB!="DRB_10minR"] ~ time[DRB!="DRB_10minR"])$coef[2],
                                                       ifelse(order == "234", lm(wf[DRB!="DRB_5minR"] ~ time[DRB!="DRB_5minR"])$coef[2],
                                                              lm(wf~time)$coef[2]))))),
         adj.r.squared = ifelse(order == "1234", summary(lm(wf~time))$adj.r.squared,
                                ifelse(order == "123", summary(lm(wf[DRB!="DRB_20minR"]~time[DRB!="DRB_20minR"]))$adj.r.squared,
                                       ifelse(order == "124", summary(lm(wf[DRB!="DRB_15minR"]~time[DRB!="DRB_15minR"]))$adj.r.squared,
                                              ifelse(order == "134", summary(lm(wf[DRB!="DRB_10minR"]~time[DRB!="DRB_10minR"]))$adj.r.squared,
                                                     ifelse(order == "234", summary(lm(wf[DRB!="DRB_5minR"]~time[DRB!="DRB_5minR"]))$adj.r.squared,
                                                            summary(lm(wf~time))$adj.r.squared))))))

### A specific case: Pol II reached the end of the gene at 10, 15 and 20min time points.
write.table(rate_df, paste0(data.dir, "elongation_rate_wave_front.txt"), quote=FALSE, sep="\t", row.names = FALSE)
```



#### test linear regression
### fit all possible linear lines with at least three experimental times
all_elongation_rate <-
  rate_df %>%
  group_by(gene_id) %>%
  mutate(fit.1234 = lm(wf~time)$coef[2],
         fit.123 = lm(wf[DRB!="DRB_20minR"] ~ time[DRB!="DRB_20minR"])$coef[2],
         fit.124 = lm(wf[DRB!="DRB_15minR"] ~ time[DRB!="DRB_15minR"])$coef[2],
         fit.134 = lm(wf[DRB!="DRB_10minR"] ~ time[DRB!="DRB_10minR"])$coef[2],
         fit.234 = lm(wf[DRB!="DRB_5minR"] ~ time[DRB!="DRB_5minR"])$coef[2]) %>%
  as.data.frame

all_elongation_rate <- all_elongation_rate[,seq(ncol(all_elongation_rate)-4, ncol(all_elongation_rate))]

all_adj.r.squared <-
  rate_df %>%
  group_by(gene_id) %>%
  mutate(fit.1234 = summary(lm(wf~time))$adj.r.squared,
         fit.123 = summary(lm(wf[DRB!="DRB_20minR"]~time[DRB!="DRB_20minR"]))$adj.r.squared,
         fit.124 = summary(lm(wf[DRB!="DRB_15minR"]~time[DRB!="DRB_15minR"]))$adj.r.squared,
         fit.134 = summary(lm(wf[DRB!="DRB_10minR"]~time[DRB!="DRB_10minR"]))$adj.r.squared,
         fit.234 = summary(lm(wf[DRB!="DRB_5minR"]~time[DRB!="DRB_5minR"]))$adj.r.squared) %>%
  as.data.frame

all_adj.r.squared <- all_adj.r.squared[,seq(ncol(all_adj.r.squared)-4, ncol(all_adj.r.squared))]
all_adj.r.squared[is.na(all_adj.r.squared)] <- 0

max_fits <- colnames(all_adj.r.squared)[(unlist(apply(all_adj.r.squared, 1, which.max)))]

### take the elongation rates with the best linear fit
all_elongation_rate %<>%
  mutate(best.fit.order = max_fits,
         gene_id = rate_df$gene_id,
         DRB = rate_df$DRB)

all_elongation_rate %<>%
  filter(DRB == "DRB_5minR") %>%
  pivot_longer(cols=-c(best.fit.order, gene_id, DRB), names_to = "fit.order", values_to = "elongation.rate") %>%
  filter(fit.order == best.fit.order) %>%
  as.data.frame

all_adj.r.squared %<>%
  mutate(best.fit.order = max_fits,
         gene_id = rate_df$gene_id,
         DRB = rate_df$DRB)

all_adj.r.squared %<>%
  filter(DRB == "DRB_5minR") %>%
  pivot_longer(cols=-c(best.fit.order, gene_id, DRB), names_to = "fit.order", values_to = "adj.r.squared") %>%
  filter(fit.order == best.fit.order) %>%
  as.data.frame

rate_df %<>%
  filter(DRB == "DRB_5minR") %>%
  left_join(all_elongation_rate %>%
              select(gene_id, best.fit.order, elongation.rate) %>%
              rename(best.elongation.rate = elongation.rate)) %>%
  left_join(all_adj.r.squared %>%
              select(gene_id, best.fit.order, adj.r.squared) %>%
              rename(best.adj.r.squared = adj.r.squared)) %>%
  dplyr::select(-c(DRB, time)) %>%
  filter(best.adj.r.squared > 0) %>%
  droplevels()

rate_df %<>%
  filter(!is.na(best.adj.r.squared), best.elongation.rate > 0)

write.table(rate_df, paste0(data.dir, "elongation_rate_wave_front.txt"), quote=FALSE, sep="\t", row.names = FALSE)