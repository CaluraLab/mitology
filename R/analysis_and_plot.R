plot_params <- function(database){
  if (database == "MitoCarta") {DBcol = c("#823cef", "#2f0967"); DBwidth = 1.1; stardist = 225
  } else if (database == "Reactome") {DBcol = c("#36ABC9", "#1d5b6b"); DBwidth = 1.3; stardist = 17
  } else if (database == "GO-CC") {DBcol = c("#ECC300", "#866f00"); DBwidth = 1.6; stardist = 14
  } else if (database == "GO-BP") {DBcol = c("#34CB80", "#1b6a43"); DBwidth = 1.3; stardist = 22}
  return(list(DBcol=DBcol, DBwidth=DBwidth, stardist=stardist))}
plot_outer_names <- function(outer_names, database, dbtree, dbdend){
  if (database == "MitoCarta") {DBoffset = 130; DBoftext = 10
    hjusts <- c(0,0,0.5,1,0,0,0)
  } else if (database == "Reactome") {DBoffset = 10; DBoftext = 2
    hjusts <- c(0,0,0,0.5,1,1,1,1,1,1,1,1,0.5,0,0,0,0,0,0,0)
  } else if (database == "GO-CC") {DBoffset = 10; DBoftext = 0.5
    hjusts <- c(0,1,0)
  } else if (database == "GO-BP") {DBoffset = 14; DBoftext = 0.5
    hjusts <- c(0,1,1,0.5,0,0,0,0,0,0,0,0)}

  if(outer_names=="leaves"){dbdend <- dbdend + ggtree::geom_tiplab(size = 1.3, offset = DBoffset)
  } else {
    Cat_nodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
    for(i in seq_along(Cat_nodes)){
      dbdend <- dbdend + ggtree::geom_cladelabel(node = Cat_nodes[i], offset = DBoffset,
                                     label = dbtree$node.label[Cat_nodes[i]-length(dbtree$tip.label)],
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



mitoanalysis_leaves <- function(dataset = NULL, database = "MitoCarta"){
  if(is.null(dataset)){stop("You must provide a dataset to perform the analysis")}
  if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
    stop("Database have to be MitoCarta, Reactome, GO-CC or GO-BP")}
  geneset <- DBgeneset(database = database)
  results <- suppressWarnings(GSVA::gsva(expr = dataset, gset.idx.list = geneset,
                                         method = "ssgsea", kcdf = "Poisson", mx.diff = F,
                                         verbose = T, ssgsea.norm = T))
  return(results)}

mitoanalysis_sections <- function(dataset = NULL, database = "MitoCarta"){
  if(is.null(dataset)){stop("You must provide a dataset to perform the analysis")}
  if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
    stop("Database have to be MitoCarta, Reactome, GO-CC or GO-BP")}
  geneset <- DBgeneset(database = database)
  dbtree <- DBtree(database = database)
  Cat_nodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
  sections <- lapply(Cat_nodes, function(x){
    clade <- ape::extract.clade(dbtree, node = x)
    unique(unname(unlist(geneset[clade$tip.label])))})
  names(sections) <- dbtree$node.label[Cat_nodes-length(dbtree$tip.label)]
  results <- suppressWarnings(GSVA::gsva(expr = dataset, gset.idx.list = sections,
                                         method = "ssgsea", kcdf = "Poisson", mx.diff = F,
                                         verbose = T, ssgsea.norm = T))
  return(results)}


plotleaves <- function(mitoanalysis_results = NULL, database = "MitoCarta",
                       samples = NULL, star = FALSE, outer_names = "sections"){

  if(is.null(mitoanalysis_results)){stop("You must provide the result of mitoanalysis_leaves to plot")}
  if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
    stop("Database have to be MitoCarta, Reactome, GO-CC or GO-BP")}
  if(outer_names != "leaves" & outer_names != "sections"){
    stop("Outer names must be leaves or sections")}

  results_score <- SummarizedExperiment::assay(mitoanalysis_results)
  dbdend <- DBdend(database = database)
  dbtree <- DBtree(database = database)

  if(!all(rownames(results_score) %in% dbtree$tip.label)){
    stop("Database does not correspond with the one used in mitoanalysis_leaves()")}
  if(is.null(samples)){samples = colnames(results_score)
  } else if(length(samples)<2){stop("Samples must be at least 2")
  } else if(!(all(samples %in% colnames(results_score)))){
    stop("All samples must be included in mitoanalysis_results")}

  plotpar <- plot_params(database = database)

  dbdend2 <- plot_outer_names(database = database, outer_names = outer_names,
                            dbtree = dbtree, dbdend = dbdend)

  g <- suppressMessages(ggtree::gheatmap(dbdend2, results_score[,colnames(results_score) %in% samples],
           colnames = T, colnames_angle = 30, font.size = 2,
           colnames_position = "top", hjust = 0, width = plotpar$DBwidth) +
    ggplot2::scale_fill_gradientn(colours = c("#701d46","#CB347F", "white", plotpar$DBcol),
                         name = "ssGSEA\nScore",
                         values = scales::rescale(c(min(results_score), 0, max(results_score))),
                         limits = c(min(results_score), max(results_score)))) +
    ggplot2::theme(legend.position = "bottom")

  return(g)
}

plotsections <- function(section_results = NULL, database = "MitoCarta", samples = NULL){

  if(is.null(section_results)){stop("You must provide the result of mitoanalysis_sections to plot")}
  if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
    stop("Database have to be MitoCarta, Reactome, GO-CC or GO-BP")}

  dbdend <- DBdend(database = database)
  dbtree <- DBtree(database = database)
  results_score <- SummarizedExperiment::assay(section_results)
  Cat_nodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
  extended_sections <- do.call(rbind, lapply(seq_along(Cat_nodes), function(x){
    clade <- ape::extract.clade(dbtree, node = Cat_nodes[x])
    results_score[rep(x, length(clade$tip.label)),]}))
  rownames(extended_sections) <- dbtree$tip.label

  # if(!all(rownames(results_score) %in% dbtree$tip.label)){
  #   stop("Database does not correspond with the one used in mitoanalysis_sections()")}
  if(is.null(samples)){samples = colnames(results_score)
  } else if(length(samples)<2){stop("Samples must be at least 2")
  } else if(!(all(samples %in% colnames(results_score)))){
    stop("All samples must be included in section_results")}

  plotpar <- plot_params(database = database)
  dbdend2 <- plot_outer_names(database = database, outer_names = "sections",
                            dbtree = dbtree, dbdend = dbdend)
  g <- suppressMessages(ggtree::gheatmap(dbdend2, extended_sections[,colnames(extended_sections) %in% samples],
                                         colnames = T, colnames_angle = 30, font.size = 2, color = NULL,
                                         colnames_position = "top", hjust = 0, width = plotpar$DBwidth) +
                          ggplot2::scale_fill_gradientn(colours=c("#701d46","#CB347F","white",plotpar$DBcol),
                                                        name = "ssGSEA\nScore",
                                                        values = scales::rescale(c(-abs(min(results_score)), 0, abs(max(results_score)))),
                                                        limits=c(-abs(min(results_score)), abs(max(results_score)))) +
                          ggplot2::theme(legend.position = "bottom"))

  return(g)
}
