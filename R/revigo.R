revigoPlots <- function(gsea_results,
                        ontol = c("BP", "MF", "CC"),
                        plotkind = c("tree", "scatter", "heat"),
                        scoretype = c("pvalue", "qvalue", "nes"),
                        plttitle = "",
                        algo = "Resnik",
                        sim_threshold = 0.7){

  res <- gseaObj@result
  res <- dplyr::filter(res, ONTOLOGY == ontol)

  # create similarity matrix
  similarity <- rrvgo::calculateSimMatrix(
    x = res$ID,
    orgdb = "org.Mm.eg.db",
    ont = ontol,
    method = algo
  )

  # set the score for each pathway using either normalized enrichment score, or normalized p or q values
  if (scoretype == "pvalue") {
    pathwayscore <- setNames(-log10(res$p.adjust), res$ID)
  } else if (scoretype == "qvalue") {
    pathwayscore <- setNames(-log10(res$qvalue), res$ID)
  } else if (scoretype == "nes") {
    pathwayscore <- setNames(res$NES, res$ID)
  }

  # reduce terms in similarity matrix
  reduced <- rrvgo::reduceSimMatrix(
    simMatrix = similarity,
    scores = pathwayscore,
    threshold = sim_threshold,
    orgdb = "org.Mm.eg.db"
  )
  if (plotkind == "tree") {
    plt <- treemap::treemap(
      dtf = reduced,
      index = c("parentTerm", "term"),
      vSize = "score",
      title = plttitle,
      palette = paletteer::paletteer_d("MoMAColors::Warhol"),
      lowerbound.cex.labels = 0,
      fontcolor.labels = c("#FFFFFFDD", "#00000080"),
      border.col = "#00000080",
      position.legend = "none",
      bg.labels = "#c5c5c5a6"
    )
  }
  if (plotkind == "scatter") {
    plt <- rrvgo::scatterPlot(similarity, reduced)
  }
  if (plotkind == "heat") {
    plt <- rrvgo::heatmapPlot(similarity, reduced)

  }

  return(plt)

}
