

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
    
    .consistencyCheck(
        database = database, nametype = nametype, 
        objectType = objectType, sections = sections)
    
    geneset <- .DBgeneset(database)
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


#' Circular heatmap on mitochondrial gene set tree.
#'
#' @description Given a matrix of scores, it returns a circular heatmap of
#' the mitochondrial gene sets (leaf of the database tree) or gene set 
#' groups (section of the database tree).
#'
#' @param data matrix or data.frame with samples in columns and mitochondrial 
#' gene sets in rows.
#' @param database character string saying the database used for the analysis.
#' Either one of "MitoCarta", "Reactome", "GO-CC" and "GO-BP".
#' @param sections logical. Either to keep the aggregated gene set categories 
#' or the specific gene sets. Default is FALSE.
#' @param samples character vector with the names of samples to be plotted.
#' Otherwise all samples are plotted.
#' @param labelNames character string that says to plot either the names of 
#' "sections" or "leaves".
#' @param ... other arguments passed on to the \code{\link[ggtree]{gheatmap}}
#' function.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @importFrom ggtree ggtree gheatmap geom_tiplab geom_cladelab
#' @importFrom scales rescale
#' @importFrom ape extract.clade Ntip
#' @importFrom magrittr %>%
#' @import ggplot2
#'
#' @examples
#' MClist <- getGeneSets()
#' n <- length(names(MClist)) * 5
#' rmatrix <- matrix(rnorm(n, 0), ncol = 5)
#' rownames(rmatrix) <- names(MClist)
#' colnames(rmatrix) <- paste0("Sample_", seq_len(5))
#' 
#'
#' @export
mitoTreeHeatmap <- function(
        data, database = "MitoCarta", sections = FALSE, 
        samples = NULL, labelNames = "sections", ...){

    .consistencyCheck(
        database = database, sections = sections, labelNames = labelNames)

    dbtree <- .DBtree(database = database)
    dbdend <- ggtree(
        tr = dbtree, layout = "fan", open.angle = 25, color = "gray60")

    if(is.null(samples)){samples <- colnames(data)
    } else if(length(samples)<2){stop("Samples must be at least 2")
    } else if(!(all(samples %in% colnames(data)))){
        stop("All samples must be included in data")}
    data <- data[,colnames(data) %in% samples]
    CatNodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
    if(sections){
        data <- do.call(rbind, lapply(seq_along(CatNodes), function(x){
            clade <- extract.clade(dbtree, node = CatNodes[x])
            data[rep(x, length(clade$tip.label)),]}))
        rownames(data) <- dbtree$tip.label }
    plotPar <- .plotParams(database = database)
    dots <- list(...)
    args <- .matchArguments(dots, list(
        p = dbdend, data = data, colnames = TRUE, 
        colnames_angle = 30, colnames_position = "top", hjust = 0, offset = 0))
    g <- do.call(gheatmap, args)
    maxScore <- max(abs(data))
    g <- g + scale_fill_gradientn(
        colours = c("#701d46","#CB347F", "white", plotPar$DBcol),
        name = "ssGSEA\nScore", limits = c(-maxScore, maxScore), 
        values = rescale(c(-maxScore, 0, maxScore)), na.value = "gray90") +
        theme(legend.position = "bottom")
    dbwidth <- (dbdend$data$x %>% range(na.rm=TRUE) %>% diff) / ncol(data)
    if(labelNames == "leaves" & !sections){
        g <- g + geom_tiplab(offset = dbwidth*(ncol(data)+1))
    } else {
        CatNames <- dbtree$node.label[CatNodes-Ntip(dbtree)]
        g <- g + geom_cladelab(
            node = CatNodes, label = CatNames, barsize = 2, 
            offset.text = (g$data$x[1]/5)/4, 
            offset = (g$data$x[1]/5)*(ncol(data)+1), 
            barcolour = plotPar$DBcol[1], angle = "auto", horizontal = TRUE) +
            theme(legend.key.size = unit(5, "mm")) }
    return(g)
}


#' Heatmap of mitochondrial gene sets.
#'
#' @description Given a matrix of scores, it returns a heatmap of
#' the mitochondrial gene sets.
#'
#' @param data matrix or data.frame with samples in columns and mitochondrial 
#' gene sets in rows.
#' @param database character string saying the database used for the analysis.
#' Either one of "MitoCarta", "Reactome", "GO-CC" and "GO-BP".
#' @param sampleAnnot character vector with samples' annotation.
#' @param splitSamples logical. If TRUE it splits samples by annotation.
#' sampleAnnot must be provided.
#' @param splitSections logical. If TRUE it splits gene sets by main section.
#' @param ... other parameters specific of the function
#' \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap-class}} object.
#'
#' @importFrom ape extract.clade
#' @importFrom ComplexHeatmap Heatmap rowAnnotation HeatmapAnnotation
#' @importFrom circlize colorRamp2
#'
#' @examples
#' MClist <- getGeneSets()
#' n <- length(names(MClist)) * 5
#' rmatrix <- matrix(rnorm(n, 0), ncol = 5)
#' rownames(rmatrix) <- names(MClist)
#' colnames(rmatrix) <- paste0("Sample_", seq_len(5))
#' mitoHeatmap(data = rmatrix, database = "MitoCarta")
#'
#' @export
mitoHeatmap <- function(
        data, database = "MitoCarta", sampleAnnot = NULL, 
        splitSamples = FALSE, splitSections = FALSE, ...){
    
    if(!(database %in% c("MitoCarta", "Reactome", "GO-CC", "GO-BP"))){
        stop("Database has to be MitoCarta, Reactome, GO-CC or GO-BP")}
    if (!is.null(sampleAnnot)) { if (length(sampleAnnot)!=ncol(data)) {
        stop("sampleAnnot length is different than samples dimension")}
    } else { if (splitSamples) { stop(
        "splitSamples can be TRUE only if sampleAnnot is provided")}}
    
    dbtree <- .DBtree(database = database)
    
    CatNodes <- dbtree$edge[,2][dbtree$edge[,1]==length(dbtree$tip.label)+1]
    CatAnnot <- unlist(lapply(seq_along(CatNodes), function(x){
        clade <- extract.clade(dbtree, node = CatNodes[x])
        rep(clade$node.label[1], length(clade$tip.label))}))
    CatAnnot <- CatAnnot[match(rownames(data), dbtree$tip.label)]
    
    plotPar <- .plotParams(database = database)
    mycol <- colorRamp2(
        c(min(data), mean(c(min(data),0)), 0, mean(c(max(data),0)), max(data)), 
        c("#701d46","#CB347F", "white", plotPar$DBcol))
    
    dots <- list(...)
    args <- .matchArguments(dots, list(
        matrix = data, name = "score", show_row_names = FALSE, col = mycol))
    
    if(splitSections){
        args$row_split <- CatAnnot
        args$row_title_rot <- 0
    } else {
        ha <- rowAnnotation(Section = CatAnnot)
        args$right_annotation <- ha }
    
    if (!is.null(sampleAnnot)) {
        if (splitSamples) {
            args$column_split <- sampleAnnot
        } else {
            hatop <- HeatmapAnnotation(sampleAnnot = sampleAnnot)
            args$top_annotation <- hatop}}
    
    g <- do.call(Heatmap, args)
    return(g)
}


#' Mitochondrial Enrichment Analysis of a gene list.
#' 
#' @description Given a vector of genes, this function will return the 
#' enrichment analysis for the mitochondrial gene sets after FDR control. For 
#' the Reactome, GO-CC and GO-BP databases it returns also the enrichment 
#' results for the corresponding original pathways.
#'
#' @param genes a vector of gene ENSEMBL id.
#' @param database character string saying the database to use for the analysis.
#' Either one of "MitoCarta", "Reactome", "GO-CC" and "GO-BP".
#'
#' @importFrom clusterProfiler enricher enrichGO
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom ReactomePA enrichPathway
#' 
#' @return enrichment analysis for the mitochondrial gene sets.
#' 
#' @export
enrichMito <- function(genes, database){
    DB_df <- getGeneSets(database = database, objectType = "dataframe")
    orM <- enricher(gene = genes, TERM2GENE = DB_df, universe = MitoGenes)
    orM <- orM@result
    
    if(database == "MitoCarta"){ return(list(MitoCarta = orM)) }
    
    if(database == "Reactome"){
        genes <- mapIds(
            x = org.Hs.eg.db, keys = genes, 
            column = "ENTREZID", keytype = "ENSEMBL")
        orO <- enrichPathway(gene = genes, maxGSSize = 1000, readable = TRUE)
        orO <- orO@result
        orO$Description <- .mitologyName(names = orO$Description, database)
    } else {
        orO <- enrichGO(
            gene = genes, keyType = "ENSEMBL", OrgDb = "org.Hs.eg.db", 
            ont = substring(database, 4, 5), pAdjustMethod = "BH", 
            readable = TRUE, pvalueCutoff = 1.1, qvalueCutoff = 1.1,
            minGSSize = 10, maxGSSize = 100000)
        orO <- orO@result
        orO$Description <- .mitologyName(names = orO$Description, "GO")
    }
    
    orO <- orO[orO$Description %in% unique(DB_df$term),]
    or <- list(orM, orO)
    names(or) <- c("mitology", "original")
    return(or)
}


#' Mitochondrial GSEA of a gene list.
#' 
#' @description Gene set enrichment analysis for the mitochondrial gene sets. 
#' For the Reactome, GO-CC and GO-BP databases it returns also the GSEA 
#' results for the corresponding original pathways.
#'
#' @param genes order ranked gene vector named by ENSEMBL id.
#' @param database character string saying the database to use for the analysis.
#' Either one of "MitoCarta", "Reactome", "GO-CC" and "GO-BP".
#'
#' @importFrom clusterProfiler GSEA gseGO
#' @importFrom AnnotationDbi mapIds
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom ReactomePA gsePathway
#' 
#' @return GSEA results for the mitochondrial gene sets.
#' 
#' @export
gseaMito <- function(genes, database){
    DB_df <- getGeneSets(database = database, objectType = "dataframe")
    orM <- GSEA(geneList = genes, TERM2GENE = DB_df)
    orM <- orM@result
    
    if(database == "MitoCarta"){ return(list(MitoCarta = orM)) }
    
    if(database == "Reactome"){
        names(genes) <- mapIds(
            x = org.Hs.eg.db, keys = names(genes), 
            column = "ENTREZID", keytype = "ENSEMBL")
        genes <- genes[!is.na(names(genes))]
        genes <- genes[!duplicated(names(genes))]
        orO <- gsePathway(geneList = genes, maxGSSize = 1000)
        orO <- orO@result
        orO$Description <- .mitologyName(names = orO$Description, database)
    } else {
        orO <- gseGO(
            geneList = genes, keyType = "ENSEMBL", OrgDb = "org.Hs.eg.db", 
            ont = substring(database, 4, 5), pAdjustMethod = "BH", 
            pvalueCutoff = 1.1, minGSSize = 10, maxGSSize = 100000)
        orO <- orO@result
        orO$Description <- .mitologyName(names = orO$Description, "GO")
    }
    orO <- orO[orO$Description %in% unique(DB_df$term),]
    or <- list(orM, orO)
    names(or) <- c("mitology", "original")
    return(or)
}


#' Circular dotplot on mitochondrial gene set tree.
#'
#' @description A circular dotplot of the mitochondrial enrichment results.
#' 
#' @param data named list of the result from enrichMito or gseaMito.
#' @param database character string saying the database to use for the analysis.
#' Either one of "MitoCarta", "Reactome", "GO-CC" and "GO-BP".
#' @param pvalCutoff pvalue cutoff to select enriched gene sets
#' @param labsize label size
#' @param max_point_size max point size
#' @param color variable used to color enriched terms, e.g. 'pvalue', 
#' 'p.adjust' or 'NES'.
#' 
#' @importFrom ggtree ggtree geom_tippoint geom_tiplab
#' @import ggplot2
#' 
#' @export
mitoTreePoint <- function(
    data, database = "MitoCarta", pvalCutoff = 0.05, 
    labsize = 3, max_point_size = 4, color = "p.adjust"){
    
    .consistencyCheck(database = database)
    plotPar <- .plotParams(database = database)
    dbtree <- .DBtree(database = database)
    g <- ggtree(tr = dbtree, layout = "fan", open.angle = 25, color = "gray60")
    
    sampleGroup <- names(data)
    data <- lapply(data, function(x) x[x$Description %in% dbtree$tip.label,])
    data <- lapply(data, function(x) x[x$p.adjust < pvalCutoff,])
    if(sum(unlist(lapply(data, function(x) nrow(x))))<1){
        stop("There are no gene sets with a p.adjust under the pvalCutoff")}
    names(data) <- sampleGroup
    
    for(i in seq_along(data)){
        g$data[, paste(color, i, sep = "_")] <- NA
        g$data[match(data[[i]]$Description, dbtree$tip.label), 
               paste(color, i, sep = "_")] <- data[[i]][,color]
        g$data[, paste("count", i, sep = "_")] <- NA
        g$data[match(data[[i]]$Description, dbtree$tip.label), 
               paste("count", i, sep = "_")] <- data[[i]][
                   ,grepl(pattern = "Count", colnames(data[[i]])) | 
                       grepl(pattern = "setSize", colnames(data[[i]]))]}
    
    xpos <- g$data$x[1] + (g$data$x[1]/5)*(seq_along(data)-1)
    
    for(i in seq_along(data)){
        g <- g + geom_tippoint(aes(
            colour = .data[[paste(color, i, sep = "_")]], 
            size = .data[[paste("count", i, sep = "_")]]), 
            x = xpos[i])}
    
    g <- g + geom_text(aes(
        x = c(xpos, rep(NA, nrow(g$data)-length(data))), 
        label = c(sampleGroup, rep(NA, nrow(g$data)-length(data))), 
        y = 0), angle = -45, hjust = 0, size = labsize)
    g$data[-match(unique(unlist(lapply(data, function(x) x$Description))), 
                  dbtree$tip.label), "label"] <- NA
    g <- g + geom_tiplab(
        offset = (g$data$x[1]/5)*(length(data)-.6), size = labsize)
    
    g <- g + scale_size_continuous(name = "", range = c(1, max_point_size)) +
        scale_colour_viridis_c(name = color) +
        theme_minimal() +
        theme(legend.position = "bottom", 
              axis.line.y = element_blank(), axis.text = element_blank(), 
              panel.grid.minor.y = element_blank(), 
              panel.grid.major.y = element_blank(), 
              panel.grid.major.x = element_line(colour = c(
                  NA, rep("grey90", length(xpos)-1), NA))) + 
        scale_x_continuous(breaks = xpos)
    
    return(g)
}