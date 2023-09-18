setwd("~/GitHub/analyzeRNAseq")
raw <- read.csv("data/counts.csv", header = TRUE, sep = ",")
rownames(raw) <- raw$Geneid
raw <- dplyr::select(raw, -Geneid)
meta <- read.csv("data/metadata.csv", header = TRUE, sep = ",")
meta$Group <- factor(meta$Group,
                         levels = c("Naive",
                                    "Tolerant",
                                    "Rejecting",
                                    "Listeria"))
rownames(meta) <- meta$SampleID
meta <- dplyr::select(meta, -SampleID)


sample_labels <- meta$Group
condition_labels <- unique(sample_labels)
condition_labels <- sort(condition_labels)
palette <- set_aes(condition_labels,
                  c("circle", "square", "triangle up", "x"),
                  c("blue", "green", "red", "orange"),
                  c("blue", "green", "red", "orange"))
anno_colors <- list(Group = c(Naive = "blue",
                                  Tolerant = "green",
                                  Rejecting = "red",
                                  Listeria = "orange"))
dseqobj <- getDESeq(raw, meta, ~Group)

rlognorm <- DESeq2::rlog(dseqobj, blind = FALSE)
rlogDF <- BiocGenerics::as.data.frame(SummarizedExperiment::assay(rlognorm))

tol_vs_naive <-  pairwiseDEGresults("Tolerant", "Naive", dseqobj)
rej_vs_tol <- pairwiseDEGresults("Rejecting", "Tolerant", dseqobj)
rej_vs_naive <- pairwiseDEGresults("Rejecting", "Naive", dseqobj)
lis_vs_tol <- pairwiseDEGresults("Listeria", "Tolerant", dseqobj)

allgenes <- tibble::rownames_to_column(rlogDF, var = "GeneID")
gmap <- createGenemap(allgenes, GMattributes = c("ensembl_gene_id", "mgi_symbol", "mgi_description"))

anno_tolvsnaive <- annotate_biomart(tol_vs_naive, gmap, c("MGI_Symbol", "MGI_Desc"))
anno_rejvsnaive <- annotate_biomart(rej_vs_naive, gmap, c("MGI_Symbol", "MGI_Desc"))
anno_rejvstol <- annotate_biomart(rej_vs_tol, gmap, c("MGI_Symbol", "MGI_Desc"))
anno_lisvstol <- annotate_biomart(lis_vs_tol, gmap, c("MGI_Symbol", "MGI_Desc"))
anno_rlogdf <- annotate_biomart(allgenes, gmap, c("MGI_Symbol", "MGI_Desc"))

