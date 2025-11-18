#!/usr/bin/env Rscript

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))
options(stringsAsFactors=FALSE)
theme_set(ggpubr::theme_pubr() + theme(legend.position="right"))

# helper variables
cols.reg <- c(Up="#e41a1c", None="grey", Down="#377eb8")
sample.days <- c("Baseline"=0, "3 DPI"=3, "4 DPI"=4, "5 DPI"=5)
genes.selected <- c("NiV"="#e41a1c", 
                    "C5"="#1f78b4", "C9"="#a6cee3", 
                    "CD209"="#33a02c", 
                    "FGB"="#6a3d9a", "FGG"="#cab2d6")
pathways.selected <- c("IFN α & β signaling"="Interferon alpha/beta signaling",
                       "APC response"="Immune response of antigen presenting cells",
                       "Macrophage activation"="Activation of macrophages",
                       "Fibrin & clotting"="Formation of Fibrin Clot (Clotting Cascade)",
                       "Complement cascade"="Complement cascade")

# helper functions
format.counts.nanostring <-  function(counts.matrix, threshold=20, nsamples=3) {
  # remove genes quantified above threshold in few samples
  keepgenes <- (counts.matrix > threshold)
  keepgenes <- rowSums(keepgenes)
  keepgenes <- (keepgenes >= nsamples)
  counts.matrix <- counts.matrix[keepgenes, ]
  
  # create DGElist and calculate CPM
  counts.matrix %>%
    DGEList() %>%
    calcNormFactors()
}
format.counts.rnaseq <- function(counts.matrix, meta.matrix, 
                          mincounts=10, minsamples=3) {
  # align count columns to metadata rows
  counts.matrix <- counts.matrix[, rownames(meta.matrix)]
  # count the number of samples (columns) per gene (row) with >= minimum counts
  x <- rowSums(counts.matrix >= mincounts)
  # keep genes (rows) with >= minimum samples above minimum counts
  x <- (x >= minsamples)
  # subset the counts matrix
  counts.matrix <- counts.matrix[x, ]
  # set the mode to integer and return
  mode(counts.matrix) <- "integer"
  return(counts.matrix)
}
get.deseq.results <- function(result.name, de.obj, psig=0.05, fsig=2) {
  # extract the comparison name
  timepoint <- result.name %>%
               str_extract("(?<=Timepoint_).+(?=_vs)") %>% 
               str_replace("\\.", " ")
  # get shrunken log2FoldChange
  rm <- results(de.obj, name=result.name, alpha=psig, cooksCutoff=FALSE)
  rm <- lfcShrink(de.obj, res=rm, coef=result.name, 
                  type="apeglm", quiet=TRUE) %>%
        as.data.frame() %>%
        rownames_to_column("Gene") %>%
        mutate(Timepoint=timepoint,
               psig=(!is.na(padj) & (padj < psig)),
               fsig=(abs(log2FoldChange) > fsig),
               Significant=(psig & fsig)) %>%
        select(-psig, -fsig)
  # format the regulation column
  rm$Regulation <- "None"
  rm[rm$Significant & (rm$log2FoldChange > 0), "Regulation"] <- "Up"
  rm[rm$Significant & (rm$log2FoldChange < 0), "Regulation"] <- "Down"
  rm$Regulation <- factor(rm$Regulation, levels=c("Up", "None", "Down"))
  # return the formatted results matrix
  return(rm)
}
format.ipa.input <- function(results.matrix) {
  results.matrix %>%
    dplyr::rename(lfc=log2FoldChange) %>%
    select(Gene, DPI, lfc, padj) %>%
    reshape2::melt(id.vars=c("Gene", "DPI")) %>% 
    mutate(variable=paste0("DPI", DPI, ".", variable)) %>%
    reshape2::dcast(Gene ~ variable, value.var="value")
}
plot.degs <- function(results.matrix, title=NULL) {
  results.matrix %>%
    filter(Significant) %>%
    group_by(DPI, Regulation) %>%
    summarise(Genes=n(),
              .groups="drop") %>%
    # add missing values
    right_join(expand.grid(Regulation=c("Up", "Down"),
                           DPI=sample.days),
               by=c("Regulation", "DPI")) %>%
    replace_na(list(Genes=0)) %>%
    # plot it
    ggplot(aes(DPI, Genes)) +
    geom_line(aes(linetype=Regulation)) +
    geom_point(aes(fill=Regulation), pch=21, size=3) +
    scale_fill_manual(values=cols.reg) +
    scale_y_continuous(limits=c(0, 725), breaks=c(0, 200, 400, 600)) +
    scale_x_continuous(breaks=sample.days, 
                       labels=names(sample.days)) +
    labs(x=NULL,
         linetype=NULL,
         fill=NULL,
         y="Total DE genes",
         title=title) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position=c(0.25, 0.8))
}
plot.genes <- function(results.matrix, genelist=genes.selected, title=NULL) {
  results.matrix %>%
    filter(Gene %in% names(genes.selected)) %>%
    select(Gene, log2FoldChange, DPI) %>%
    # add baseline
    rbind(data.frame(Gene=names(genes.selected),
                     log2FoldChange=0,
                     DPI=0)) %>%
    ggplot(aes(DPI, log2FoldChange, col=Gene)) +
    geom_hline(yintercept=2, linetype=2, col="grey60") +
    geom_line() +
    geom_point() +
    scale_color_manual(values=genes.selected) +
    scale_x_continuous(breaks=sample.days) +
    scale_y_continuous(limits=c(NA, 18),
                       breaks=c(0, 2, 5, 10, 15)) +
    labs(x=NULL,
         col=NULL,
         y="Fold change (log2)",
         title=title) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position="none")
}
plot.pathways <- function(filename, pathways=pathways.selected, title=NULL) {
  filename %>%
    read.csv() %>%
    select(-Type) %>%
    reshape2::melt(id.vars="Name",
                   variable.name="DPI",
                   value.name="zscore") %>%
    # format name and DPI
    filter(Name %in% pathways) %>%
    mutate(Name=factor(Name, 
                       levels=pathways, 
                       labels=names(pathways)),
           DPI=as.integer(str_extract(DPI, "[0-9]+"))) %>%
    # fill in baseline
    rbind(data.frame(Name=names(pathways),
                     DPI=0, 
                     zscore=0)) %>%
    # fill in missing values and NAs with zeros (NA or missing = not sig.)
    reshape2::dcast(Name ~ DPI, value.var="zscore", fill=0) %>%
    reshape2::melt(id.vars="Name",
                   variable.name="DPI",
                   value.name="zscore") %>% 
    mutate(DPI=as.integer(as.character(DPI))) %>% 
    ggplot(aes(DPI, zscore, col=Name)) +
    geom_line(aes(group=Name)) +
    geom_point() +
    scale_color_brewer(palette="Set2") +
    scale_x_continuous(breaks=sample.days) +
    scale_y_continuous(limits=c(NA, 5.1)) +
    labs(x=NULL,
         col=NULL,
         y="Pathway z-score",
         title=title) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          legend.position="none")
} 

# metadata: remove samples that didn't pass QC
meta <- read.csv("samplesheet.csv", na.strings="") %>%
        filter(QC.pass) %>%
        mutate(Timepoint=factor(Timepoint, levels=names(sample.days)))
row.names(meta) <- meta$ID

## Nanostring ------------------------------------------------------------------
# subset metadata table
m <- filter(meta, Assay=="nanoString")

# read and format thresholded matrix
cmat <- read.csv("data/nanostring-thresholded.csv") %>%
        # remove positive and negative controls
        filter(!str_detect(Gene, "^NEG_|^POS_")) %>%
        column_to_rownames("Gene") %>%
        as.matrix()
# align rows and columns
x <- intersect(rownames(meta), colnames(cmat))
m <- m[x, ]
cmat <- cmat[, x]
rm(x)

# define a model matrix with "Experiment" as batch and "Timepoint" as group
modmat <- m %>% # formatting factor for easy column names
          mutate(Timepoint=factor(Timepoint,
                                  levels=names(sample.days),
                                  labels=c("Baseline", "DPI3", "DPI4", "DPI5")))
modmat <- model.matrix(~0 + Timepoint + Experiment, data=modmat)

# get normalized counts and integrate the model
rmat <- cmat %>%
        format.counts.nanostring() %>%
        voom(modmat) %>%
        lmFit(modmat)
# define contrasts and add to the model
cont <- makeContrasts(TimepointDPI3 - TimepointBaseline, 
                      TimepointDPI4 - TimepointBaseline,
                      TimepointDPI5 - TimepointBaseline,
                      levels=colnames(modmat))
rmat <- contrasts.fit(rmat, cont) %>%
        eBayes(robust=TRUE)

# get significance and regulation
sigmat <- rmat %>%
          decideTests(lfc=1) %>%
          as.data.frame() %>%
          rownames_to_column("Gene") %>%
          reshape2::melt(id.vars="Gene",
                         variable.name="Timepoint",
                         value.name="Regulation") %>%
          mutate(Significant=(Regulation != 0),
                 Regulation=factor(Regulation, 
                                   levels=c(1, 0, -1), 
                                   labels=names(cols.reg)))

# get p-values and log2 fold changes for all genes and comparisons
rmat <- colnames(cont) %>%
        lapply(function(i) {
          topTable(rmat, coef=i, number=Inf) %>%
            rownames_to_column("Gene") %>%
            mutate(Timepoint=i)
        }) %>%
        do.call(rbind, .) %>%
        # add significance
        full_join(sigmat, by=c("Gene", "Timepoint")) %>% 
        # rename some columns for easy recall
        dplyr::rename(lfc=logFC, p=P.Value, padj=adj.P.Val) %>%
        # format time point, add DPI, and add significance columns
        mutate(Timepoint=str_extract(Timepoint, "(?<=^Timepoint)[^ ]+"),
               Timepoint=factor(Timepoint, 
                                levels=c("Baseline", "DPI3", "DPI4", "DPI5"),
                                labels=names(sample.days)),
               DPI=factor(Timepoint, 
                          levels=names(sample.days), 
                          labels=sample.days),
               DPI=as.integer(as.character(DPI))) %>%
        # only keep the columns we need
        select(Gene, lfc, p, padj, Timepoint, DPI, Significant, Regulation)
write.csv(rmat, "analysis/results-blood.csv", row.names=FALSE)

# clean up 
rm(modmat, cont, sigmat)

# plot volcano with the top DE genes labeled
topgenes <- rmat %>%
            filter(Significant) %>%
            group_by(Timepoint) %>%
            top_n(n=15, wt=abs(lfc)) %>%
            ungroup()
# all days
rmat %>%
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(size=Regulation, col=Regulation), alpha=0.5) +
  ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                           size=2, force=10, max.time=10, max.overlaps=20) +
  scale_size_manual(values=c(1, 0.5, 1)) +
  scale_color_manual(values=cols.reg) +
  facet_wrap(~Timepoint, nrow=1) +
  labs(x="Fold change (log2)",
       y="FDR-adjusted p-value (-log10)") +
  theme(legend.position="right")
ggsave("analysis/volcano-blood-all.png", 
       units="in", width=6.5, height=3)
# 3 DPI only
topgenes <- rmat %>%
            filter(Significant, DPI==3) %>%
            top_n(n=8, wt=abs(lfc)) 
rmat %>%
  filter(DPI==3) %>%
  ggplot(aes(lfc, -log10(padj))) +
  geom_point(aes(size=Regulation, col=Regulation)) +
  scale_color_manual(values=cols.reg) +
  scale_size_manual(values=c(1, 0.5, 1)) +
  ggrepel::geom_text_repel(data=topgenes, aes(label=Gene), 
                           force=10, max.time=10, max.overlaps=20) +  
  labs(x="Fold change (log2)",
       y="FDR-adj. p-value (-log10)",
       title="Blood, 3 DPI") +
  theme(legend.position="none")
ggsave("analysis/volcano-blood-3dpi.png",
       units="in", width=2.8, height=3)

# clean up
rm(topgenes, m, rmat)

## RNA-seq setup ---------------------------------------------------------------
# subset metadata
meta <- filter(meta, Assay=="RNA-seq")

# load the counts matrix
cmat <- read.delim("data/counts.tsv", 
                   sep="\t") %>%
        # remove rRNA, etc.
        filter(!str_detect(Gene, "_rRNA"),
               Gene!="Metazoa_SRP",
               Gene!="Y_RNA") %>%
        # collapse to one row per gene name
        reshape2::melt(id.vars="Gene") %>%
        group_by(Gene, variable) %>%
        summarise(value=sum(value),
                  .groups="drop") %>%
        reshape2::dcast(Gene ~ variable, value.var="value") %>%
        # format as a matrix
        column_to_rownames("Gene") %>%
        as.matrix()

# collapse viral genes to "NiV"
vgenes <- c("N", "P/V/W", "C", "M", "F", "G", "L")
# get all host genes
hgenes <- cmat[!(rownames(cmat) %in% vgenes), ]
# get the sum of all quantified viral reads for each sample
vgenes <- colSums(cmat[vgenes, ])
# format viral reads as a 1xN matrix
vgenes <- matrix(vgenes, nrow=1, dimnames=list("NiV", names(vgenes)))
cmat <- rbind(hgenes, vgenes)
rm(vgenes, hgenes)

# align rows and columns
x <- intersect(rownames(meta), colnames(cmat))
meta <- meta[x, ]
cmat <- cmat[, x]
rm(x)

## RNA-seq tonsil --------------------------------------------------------------
# subset to tonsil and run DESeq2
m <- filter(meta, Tissue=="Tonsil")
d <- format.counts.rnaseq(cmat, m)
d <- DESeqDataSetFromMatrix(d, m, ~Timepoint) %>%
     DESeq(quiet=TRUE)

# get a results matrix with shrunken log2FC
r <- resultsNames(d)[-1] %>%
     lapply(get.deseq.results, de.obj=d) %>%
     do.call(rbind, .) %>%
     # get DPI column and format the "Timepoint" column
     mutate(DPI=as.integer(str_extract(Timepoint, "^[0-9]+")),
            Timepoint=factor(Timepoint, levels=names(sample.days)))
write.csv(r, "analysis/results-tonsil.csv", row.names=FALSE)

# save IPA input
r %>%
  format.ipa.input() %>% 
  write.csv("analysis/ipa-input-tonsil.csv", row.names=FALSE)

# total DE genes
plot.degs(r, title="Tonsil DE genes")
ggsave("analysis/degenes-tonsil.png", units="in", width=2.8, height=3)

# selected genes
plot.genes(r, title="Tonsil marker genes")
ggsave("analysis/markers-tonsil.png", units="in", width=2.8, height=3)

# load IPA and plot selected pathways
plot.pathways("analysis/ipa-output-tonsil.csv", title="Tonsil pathways")
ggsave("analysis/pathways-tonsil.png", units="in", width=2.8, height=3)

# clean up
rm(d, m, r)

## RNA-seq lung ----------------------------------------------------------------
# subset to lung and run DESeq2
m <- filter(meta, Tissue=="Lung")
d <- format.counts.rnaseq(cmat, m)
d <- DESeqDataSetFromMatrix(d, m, ~Timepoint) %>%
     DESeq(quiet=TRUE)

# get a results matrix with shrunken log2FC
r <- resultsNames(d)[-1] %>%
     lapply(get.deseq.results, de.obj=d) %>%
     do.call(rbind, .) %>%
     # get DPI column and format the "Timepoint" column
     mutate(DPI=as.integer(str_extract(Timepoint, "^[0-9]+")),
            Timepoint=factor(Timepoint, levels=names(sample.days)))
write.csv(r, "analysis/results-lung.csv", row.names=FALSE)

# save IPA input
r %>%
  format.ipa.input() %>% 
  write.csv("analysis/ipa-input-lung.csv", row.names=FALSE)

# total DE genes
plot.degs(r, title="Lung DE genes")
ggsave("analysis/degenes-lung.png", units="in", width=2.8, height=3)

# selected genes
plot.genes(r, title="Lung marker genes")
ggsave("analysis/markers-lung.png", units="in", width=2.8, height=3)

# load IPA and plot selected pathways
plot.pathways("analysis/ipa-output-lung.csv", title="Lung pathways")
ggsave("analysis/pathways-lung.png", units="in", width=2.8, height=3)
  
# clean up
rm(d, m, r)

## fin -------------------------------------------------------------------------
sessionInfo()
