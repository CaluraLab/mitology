plotParams <- function(database){
  if (database == "MitoCarta") {DBcol = c("#823cef", "#2f0967"); DBwidth = 1.1; stardist = 225
  } else if (database == "Reactome") {DBcol = c("#36ABC9", "#1d5b6b"); DBwidth = 1.3; stardist = 17
  } else if (database == "GO-CC") {DBcol = c("#ECC300", "#866f00"); DBwidth = 1.6; stardist = 14
  } else if (database == "GO-BP") {DBcol = c("#34CB80", "#1b6a43"); DBwidth = 1.3; stardist = 22}
  return(list(DBcol = DBcol, DBwidth = DBwidth, stardist = stardist))}

plotOuterNames <- function(outerNames, database, dbtree, dbdend){
  if (database == "MitoCarta") {DBoffset = 130; DBoftext = 10
    hjusts <- c(0,0,0.5,1,0,0,0)
  } else if (database == "Reactome") {DBoffset = 10; DBoftext = 2
    hjusts <- c(0,0,0,0.5,1,1,1,1,1,1,1,1,0.5,0,0,0,0,0,0,0)
  } else if (database == "GO-CC") {DBoffset = 10; DBoftext = 0.5
    hjusts <- c(0,1,0)
  } else if (database == "GO-BP") {DBoffset = 14; DBoftext = 0.5
    hjusts <- c(0,1,1,0.5,0,0,0,0,0,0,0,0)}

  if(outerNames=="leaves"){dbdend <- dbdend + ggtree::geom_tiplab(size = 1.3, offset = DBoffset)
  } else {
    CatNodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
    for(i in seq_along(CatNodes)){
      dbdend <- dbdend + ggtree::geom_cladelabel(node = CatNodes[i], offset = DBoffset,
                                     label = dbtree$node.label[CatNodes[i]-length(dbtree$tip.label)],
                                     hjust = hjusts[i], offset.text = DBoftext,
                                     fontsize = 3, barsize = 0)}}
  return(dbdend)}

DBgeneset <- function(database){
  if (database == "MitoCarta") {geneset = MCgs} else if (database == "Reactome") {geneset = RTgs
  } else if (database == "GO-CC") {geneset = CCgs} else if (database == "GO-BP") {geneset = BPgs}
  return(geneset)}

DBtree <- function(database){
  if (database == "MitoCarta") {dbtree = MCtree} else if (database == "Reactome") {dbtree = RTtree
  } else if (database == "GO-CC") {dbtree = CCtree} else if (database == "GO-BP") {dbtree = BPtree}
  return(dbtree)}

DBdend <- function(database){
  if (database == "MitoCarta") {dbdend = MCdend} else if (database == "Reactome") {dbdend = RTdend
  } else if (database == "GO-CC") {dbdend = CCdend} else if (database == "GO-BP") {dbdend = BPdend}
  return(dbdend)}


#' Perform ssGSEA analysis by leaves.
#'
#' Given a dataset, it returns the ssGSEA score for each leaves and each sample.
#'
#' @param dataset Gene expression data which can be given either as a SummarizedExperiment or ExpressionSet object, or as a matrix of expression values where rows correspond to genes and columns correspond to samples.
#' @param database One of MitoCarta, Reactome, GO-CC and GO-BP. Default is MitoCarta.
#'
#' @return NULL
#'
#' @importFrom GSVA gsva
#'
#' @export
mitoAnalysisLeaves <- function(dataset = NULL, database = "MitoCarta"){

  if(is.null(dataset)){stop("You must provide a dataset to perform the analysis")}
  if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
    stop("Database have to be MitoCarta, Reactome, GO-CC or GO-BP")}

  geneset <- DBgeneset(database = database)
  results <- suppressWarnings(gsva(expr = dataset, gset.idx.list = geneset,
                                   method = "ssgsea", kcdf = "Poisson", mx.diff = F,
                                   verbose = T, ssgsea.norm = T))
  return(results)}


#' Perform ssGSEA analysis by sections.
#'
#' Given a dataset, it returns the ssGSEA score for each section and each sample.
#'
#' @inheritParams mitoAnalysisLeaves
#'
#' @return NULL
#'
#' @importFrom ape extract.clade
#' @importFrom GSVA gsva
#'
#' @export
mitoAnalysisSections <- function(dataset = NULL, database = "MitoCarta"){

  if(is.null(dataset)){stop("You must provide a dataset to perform the analysis")}
  if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
    stop("Database have to be MitoCarta, Reactome, GO-CC or GO-BP")}

  geneset <- DBgeneset(database = database)
  dbtree <- DBtree(database = database)

  CatNodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
  sections <- lapply(CatNodes, function(x){
    clade <- extract.clade(dbtree, node = x)
    unique(unname(unlist(geneset[clade$tip.label])))})
  names(sections) <- dbtree$node.label[CatNodes-length(dbtree$tip.label)]

  results <- suppressWarnings(gsva(expr = dataset, gset.idx.list = sections,
                                   method = "ssgsea", kcdf = "Poisson", mx.diff = F,
                                   verbose = T, ssgsea.norm = T))
  return(results)}


#' Plot a circular heatmap of leaf ssGSEA scores on the dendrogram.
#'
#' @param leafResults Result of the mitoAnalysisLeaves function.
#' @param database One of MitoCarta, Reactome, GO-CC and GO-BP. Default is MitoCarta.
#' @param samples A vector with the names of samples to be plotted. If NULL all samples are plotted.
#' @param outerNames Whether to plot section names or leaves names. Default is sections.
#'
#' @return NULL
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom ggtree gheatmap
#' @importFrom ggplot2 scale_fill_gradientn theme
#' @importFrom scales rescale
#'
#' @export
mitoPlotLeaves <- function(leafResults = NULL, database = "MitoCarta",
                       samples = NULL, outerNames = "sections"){

  if(is.null(leafResults)){stop("You must provide the result of mitoAnalysisLeaves to plot")}
  if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
    stop("Database have to be MitoCarta, Reactome, GO-CC or GO-BP")}
  if(outerNames != "leaves" & outerNames != "sections"){
    stop("Outer names must be leaves or sections")}

  resultsScore <- assay(leafResults)
  dbdend <- DBdend(database = database)
  dbtree <- DBtree(database = database)

  if(!all(rownames(resultsScore) %in% dbtree$tip.label)){
    stop("Database does not correspond with the one used in mitoAnalysisLeaves()")}
  if(is.null(samples)){samples = colnames(resultsScore)
  } else if(length(samples)<2){stop("Samples must be at least 2")
  } else if(!(all(samples %in% colnames(resultsScore)))){
    stop("All samples must be included in leafResults")}

  plotpar <- plotParams(database = database)

  dbdend2 <- plotOuterNames(database = database, outerNames = outerNames,
                              dbtree = dbtree, dbdend = dbdend)

  g <- suppressMessages(gheatmap(dbdend2, resultsScore[,colnames(resultsScore) %in% samples],
                                 colnames = T, colnames_angle = 30, font.size = 2,
                                 colnames_position = "top", hjust = 0, width = plotpar$DBwidth) +
                          scale_fill_gradientn(colours = c("#701d46","#CB347F", "white", plotpar$DBcol),
                                               name = "ssGSEA\nScore",
                                               values = rescale(c(min(resultsScore), 0, max(resultsScore))),
                                               limits = c(min(resultsScore), max(resultsScore)))) +
    theme(legend.position = "bottom")

  return(g)
}


#' Plot a circular heatmap of section ssGSEA scores on the dendrogram.
#'
#' @param sectionResults Result of the mitoAnalysisSections function.
#' @param database One of MitoCarta, Reactome, GO-CC and GO-BP. Default is MitoCarta.
#' @param samples A vector with the names of samples to be plotted. If NULL all samples are plotted.
#'
#' @return NULL
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom ggtree gheatmap
#' @importFrom ggplot2 scale_fill_gradientn theme
#' @importFrom scales rescale
#' @importFrom ape extract.clade
#'
#' @export
mitoPlotSections <- function(sectionResults = NULL, database = "MitoCarta", samples = NULL){

  if(is.null(sectionResults)){stop("You must provide the result of mitoAnalysisSections to plot")}
  if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
    stop("Database have to be MitoCarta, Reactome, GO-CC or GO-BP")}

  dbdend <- DBdend(database = database)
  dbtree <- DBtree(database = database)
  resultsScore <- assay(sectionResults)
  CatNodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
  extSections <- do.call(rbind, lapply(seq_along(CatNodes), function(x){
    clade <- extract.clade(dbtree, node = CatNodes[x])
    resultsScore[rep(x, length(clade$tip.label)),]}))
  rownames(extSections) <- dbtree$tip.label

  # if(!all(rownames(resultsScore) %in% dbtree$tip.label)){
  #   stop("Database does not correspond with the one used in mitoAnalysisSections()")}
  if(is.null(samples)){samples = colnames(resultsScore)
  } else if(length(samples)<2){stop("Samples must be at least 2")
  } else if(!(all(samples %in% colnames(resultsScore)))){
    stop("All samples must be included in sectionResults")}

  plotpar <- plotParams(database = database)
  dbdend2 <- plotOuterNames(database = database, outerNames = "sections",
                            dbtree = dbtree, dbdend = dbdend)
  g <- suppressMessages(gheatmap(dbdend2, extSections[,colnames(extSections) %in% samples],
                                 colnames = T, colnames_angle = 30, font.size = 2, color = NULL,
                                 colnames_position = "top", hjust = 0, width = plotpar$DBwidth) +
                          scale_fill_gradientn(colours=c("#701d46","#CB347F","white",plotpar$DBcol),
                                               name = "ssGSEA\nScore",
                                               values = rescale(c(-abs(min(resultsScore)), 0, abs(max(resultsScore)))),
                                               limits=c(-abs(min(resultsScore)), abs(max(resultsScore)))) +
                          theme(legend.position = "bottom"))

  return(g)
}
