

#' Get the mitochondrial gene sets
#'
#' @description It returns the mitochondrial gene sets (in form of list or 
#' data frame) of the four possible databases: "MitoCarta", "Reactome", 
#' "GO-CC" and "GO-BP".
#'
#' @param database character string saying the database to use for the analysis.
#' Either one of "MitoCarta", "Reactome", "GO-CC" and "GO-BP".
#' @param nametype character string saying the type of gene name ID. 
#' Either one of "SYMBOL", "ENTREZID" or "ENSEMBL".
#' @param objectType character string saying the type of needed object. 
#' Either one of "list" or "dataframe".
#' @param sections logical. Either to keep the aggregated gene set categories 
#' or the specific gene sets. Default is FALSE.
#'
#' @importFrom ape extract.clade
#' 
#' @return the mitochondrial gene sets.
#'
#' @examples
#' MClist <- getGeneSets()
#'
#' @export
getGeneSets <- function(
        database = "MitoCarta", nametype = "ENSEMBL", 
        objectType = "list", sections = FALSE){
    
    .consistencyCheck(database, nametype, objectType, sections)
    
    geneset <- .DBgeneset(database)
    # geneset <- geneset[[nametype]]
    geneset <- .geneIDtrans(nametype, geneset, database)
    
    if(sections){
        dbtree <- .DBtree(database)
        CatNodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
        geneset <- lapply(CatNodes, function(x){
            clade <- extract.clade(dbtree, node = x)
            unique(unname(unlist(geneset[clade$tip.label])))})
        names(geneset) <- dbtree$node.label[CatNodes-length(dbtree$tip.label)]}
    
    if(objectType == "dataframe"){
        geneset <- do.call(rbind, lapply(seq_along(geneset), function(i){
            data.frame(
                term = rep(names(geneset)[i], length(geneset[[i]])),
                gene = geneset[[i]]) }))}
    
    return(geneset)
}


#' Circular heatmap of mitochondrial gene set scores.
#'
#' @description Given the enrichment output, it returns a circular heatmap of 
#' the ssGSEA scores for each specific mitochondrial gene set (leaf of the 
#' database tree) or gene set group (section of the database tree).
#'
#' @param leafResults output of the enrichment.
#' @param database character string saying the database used for the analysis.
#' Either one of "MitoCarta", "Reactome", "GO-CC" and "GO-BP".
#' @param sections logical. Either to keep the aggregated gene set categories 
#' or the specific gene sets. Default is FALSE.
#' @param samples character vector with the names of samples to be plotted.
#' Otherwise all samples are plotted.
#' @param labelNames character string that says to plot the names of "sections"
#' or "leaves".
#' @param ... other arguments passed on to the \code{\link[ggtree]{gheatmap}}
#' function.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom ggtree ggtree gheatmap geom_tiplab geom_cladelab
#' @importFrom scales rescale
#' @importFrom ape extract.clade Ntip
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @examples
#' MClist <- getGeneSets()
#'
#' @export
mitoPlotLeaves <- function(
        leafResults, database = "MitoCarta", sections = FALSE, 
        samples = NULL, labelNames = "sections", ...){

    if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
        stop("Database has to be MitoCarta, Reactome, GO-CC or GO-BP")}
    if(labelNames != "leaves" & labelNames != "sections"){
        stop("Label names must be leaves or sections")}

    dbtree <- .DBtree(database = database)
    dbdend <- ggtree(
        tr = dbtree, layout = "fan", open.angle = 25, color = "gray60")

    # if(!all(rownames(leafResults) %in% dbtree$tip.label)){stop(
    #     "Database does not correspond with the results")}
    if(is.null(samples)){samples <- colnames(leafResults)
    } else if(length(samples)<2){stop("Samples must be at least 2")
    } else if(!(all(samples %in% colnames(leafResults)))){
        stop("All samples must be included in leafResults")}

    resultsScore <- leafResults[,colnames(leafResults) %in% samples]

    CatNodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
    
    if(sections){
        resultsScore <- do.call(rbind, lapply(seq_along(CatNodes), function(x){
            clade <- extract.clade(dbtree, node = CatNodes[x])
            resultsScore[rep(x, length(clade$tip.label)),]}))
        rownames(resultsScore) <- dbtree$tip.label }
    
    plotPar <- .plotParams(database = database)

    dots <- list(...)
    args <- .matchArguments(dots, list(
        p = dbdend, data = resultsScore, colnames = TRUE, 
        colnames_angle = 30, colnames_position = "top", hjust = 0, offset = 0))
    g <- do.call(gheatmap, args)
    maxScore <- max(abs(resultsScore))
    g <- g + scale_fill_gradientn(
        colours = c("#701d46","#CB347F", "white", plotPar$DBcol),
        name = "ssGSEA\nScore",
        values = rescale(c(-maxScore, 0, maxScore)),
        limits = c(-maxScore, maxScore), 
        na.value = "gray90") +
        theme(legend.position = "bottom")

    dbwidth <- (
        dbdend$data$x %>% range(na.rm=TRUE) %>% diff) / ncol(resultsScore)

    if(labelNames == "leaves" & !sections){
        g <- g + geom_tiplab(offset = dbwidth*(ncol(resultsScore)+1))
    } else {
        CatNames <- dbtree$node.label[CatNodes-Ntip(dbtree)]
        g <- g + geom_cladelab(
            node = CatNodes, label = CatNames, barsize = 2, 
            offset.text = dbwidth*2, offset = dbwidth*(ncol(resultsScore)+1), 
            barcolour = plotPar$DBcol[1], angle = "auto", horizontal = TRUE) +
            theme(legend.key.size = unit(5, "mm")) }

    return(g)
}

