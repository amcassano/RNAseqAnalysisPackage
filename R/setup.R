# setwd("~/GitHub/analyzeRNAseq")
# raw <- read.csv("data/counts.csv", header = TRUE, sep = ",")
# rownames(raw) <- raw$Geneid
# raw <- dplyr::select(raw, -Geneid)
#  meta <- read.csv("data/metadata.csv", header = TRUE, sep = ",")
#  meta$Group <- factor(meta$Group,
#                          levels = c("Naive",
#                                     "Tolerant",
#                                     "Rejecting"))
# rownames(meta) <- meta$SampleID
# meta <- dplyr::select(meta, -SampleID)
#
#
# sample_labels <- meta$Group
# condition_labels <- unique(sample_labels)
# condition_labels <- sort(condition_labels)
# palette <- set_aes(condition_labels,
#                   c("circle", "square", "triangle up"),
#                   c("blue", "green", "red"),
#                   c("blue", "green", "red"))
# anno_colors <- list(Group = c(Naive = "blue",
#                                   Tolerant = "green",
#                                   Rejecting = "red"))
# dseqobj <- getDESeq(raw, meta, ~Group)
#
# rlognorm <- DESeq2::rlog(dseqobj, blind = FALSE)
# rlogDF <- BiocGenerics::as.data.frame(SummarizedExperiment::assay(rlognorm))
#
# tol_vs_naive <-  pairwiseDEGresults("Tolerant", "Naive", dseqobj)
# rej_vs_tol <- pairwiseDEGresults("Rejecting", "Tolerant", dseqobj)
# rej_vs_naive <- pairwiseDEGresults("Rejecting", "Naive", dseqobj)
#
# allgenes <- tibble::rownames_to_column(rlogDF, var = "GeneID")
# allgenes <- select(allgenes, GeneID)
# gmap <- createGenemap(allgenes, GMattributes = c("ensembl_gene_id", "mgi_symbol", "mgi_description"))
#
# anno_tolvsnaive <- annotate_biomart(tol_vs_naive, gmap, c("MGI_Symbol", "MGI_Desc"))
# anno_rejvsnaive <- annotate_biomart(rej_vs_naive, gmap, c("MGI_Symbol", "MGI_Desc"))
# anno_rejvstol <- annotate_biomart(rej_vs_tol, gmap, c("MGI_Symbol", "MGI_Desc"))
# anno <- rownames_to_column(rlogDF, var = "GeneID")
# anno_rlogdf <- annotate_biomart(anno, gmap, c("MGI_Symbol", "MGI_Desc"))
#
# ranked_tolvsnaive <- signif_deg(anno_tolvsnaive)
# # export_genelist(ranked_tolvsnaive, "tolvsnaive")
# ranked_rejvsnaive <- signif_deg(anno_rejvsnaive)
#
