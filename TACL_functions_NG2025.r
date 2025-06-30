#TACl functions

norm4C <- function(readsGR, nReads = 1e6, nTop = 2, wSize = 21) {
  readsGR$normReads <- 0
  sumTop <- sum(-sort(-readsGR$reads)[1:nTop])
  wNorm <- nReads / (sum(readsGR$reads) - sumTop)
  readsGR$normReads <- wNorm * readsGR$reads
  readsGR$norm4C <- runmean(x = readsGR$normReads, k = wSize, endrule = "mean")
  return(readsGR)
}

local4C <- function(rds, ROI, normZoom = 2e7) {
  zoom4C <- resize(ROI, width = normZoom, fix = "center")
  GR1 <- subsetByOverlaps(rds$reads, zoom4C)
  GR1 <- subset(GR1, type == "non_blind")
  return(GR1)
}

normalize4C <- function(list4C, listNames = c("on_2", "off_2", "mch_2"), ROI, normZoom = 2e7, wSize = 21) {
  # Filter the list for specified layer names
  list4C_OI <- list4C[listNames]

  # Zoom around the ROI for each layer
  zoom4C <- resize(ROI, width = normZoom, fix = "center")
  list4C_OI <- lapply(list4C_OI, function(gr) {
    overlapped <- gr[overlapsAny(gr, zoom4C)]
  })

  # Filter each GRanges object to keep only ranges with common fe_ids
  feIdLists <- lapply(list4C_OI, function(x) unique(x$fe_id))
  commonFeIds <- Reduce(intersect, feIdLists)
  list4C_OI <- lapply(list4C_OI, function(gr) {
    gr[gr$fe_id %in% commonFeIds]
  })

  # Normalize each GRanges object locally around each TetO site
  results <- list()

  # Loop over each TetO site
  for (i in seq_along(zoom4C)) {
    message(paste0("Normalizing TetO site ", i, " of ", length(zoom4C)))
    results[[i]] <- list()

    for (layerName in names(list4C_OI)) {
      message(paste0("  Normalizing layer ", layerName))
      layerGR <- list4C_OI[[layerName]]
      overlaps <- subsetByOverlaps(layerGR, zoom4C[i])
      normalized <- norm4C(overlaps, nReads = 1e6, nTop = 2, wSize = wSize)
      results[[i]][[layerName]] <- normalized
    }
  }

  # Combine the normalized results into a single GRanges object for each TetO site
  combinedResults <- list()

  # Loop over each TetO site
  for (i in seq_along(results)) {
    tetOData <- results[[i]]
    grList <- list()

    # Extract and rename the desired columns for each layer
    for (layerName in names(tetOData)) {
      gr <- tetOData[[layerName]]

      # Select only the 'pos', 'reads', and 'norm4C' columns
      mcols(gr) <- mcols(gr)[, c("pos", "reads", "norm4C")]

      # Store the processed GRanges object in grList
      grList[[layerName]] <- gr
    }

    # Initialize combinedGR using the genomic ranges from the first item in grList
    combinedGR <- grList[[1]][, seq_len(0), drop = FALSE]
    mcols(combinedGR)$pos <- mcols(grList[[1]])$pos

    # Merge the metadata columns from all layers into combinedGR
    for (layerName in names(grList)) {
      gr <- grList[[layerName]]
      newReadsColName <- paste0("reads_", layerName)
      newNorm4CColName <- paste0("norm4C_", layerName)

      # Add the renamed columns to combinedGR's metadata
      mcols(combinedGR)[[newReadsColName]] <- mcols(gr)$reads
      mcols(combinedGR)[[newNorm4CColName]] <- mcols(gr)$norm4C
    }

    # Store the combined GRanges object for this TetO site
    combinedResults[[i]] <- combinedGR
  }

  return(combinedResults)
}

addTetOdistance<-function(peaks,C17TetO){
  distHits<-distanceToNearest(peaks, C17TetO)
  peaks$distance<-1e9 #placeholder
  peaks[queryHits(distHits)]$distance<-mcols(distHits)$distance
  return(peaks)
}

getMeanBG_V3 <- function(region, bw, signalValueScore = 0, distanceCutoff = 3e6, domain = domain_HMM_GR, binWidth = 10) {
    if (signalValueScore == 0) {
        message("No filtering for signalValue")
        region <- region[(region$distance > distanceCutoff | is.na(region$distance))]
        region <- subsetByOverlaps(region, domain, invert = TRUE)
    } else {
        region <- region[region$signalValue > signalValueScore & (region$distance > distanceCutoff | is.na(region$distance))]
        region <- subsetByOverlaps(region, domain, invert = TRUE)
    }

    region <- resize(x = region, width = binWidth, fix = "center")
    bw_score <- import(bw, selection = BigWigSelection(region))
    gr.data.cov <- GenomicRanges::coverage(bw_score, weight = "score")
    seqlevels(region) <- names(gr.data.cov)
    region <- binnedAverage(region, gr.data.cov, "bwscore")
    meanBG <- round(x = mean(region$bwscore), digits = 2)

    return(meanBG)
}

getOuterBinsBG_V3 <- function(region, bw, signalValueScore = 0, distanceCutoff = 3e6, domain, zoomWidth = 5000, flankWidth = 1000) {
    # region=peaks
    # bw=paste0(prefix,ChIPdf$bigwig[a])
    # signalValueScore=ChIPdf$signalValue[a]
    # distanceCutoff=3e6

    # start_time <- Sys.time()

    if (signalValueScore == 0) {
        message("No filtering for signalValue")
        region <- region[(region$distance > distanceCutoff | is.na(region$distance))]
        region <- subsetByOverlaps(region, domain, invert = TRUE)
    } else {
        region <- region[region$signalValue > signalValueScore & (region$distance > distanceCutoff | is.na(region$distance))]
        region <- subsetByOverlaps(region, domain, invert = TRUE)
    }



    region <- resize(x = region, width = zoomWidth, fix = "center")
    region_Up <- resize(x = region, width = flankWidth, fix = "start")
    region_Down <- resize(x = region, width = flankWidth, fix = "end")
    region <- c(region_Up, region_Down)
    bw_score <- import(bw, selection = BigWigSelection(region))
    gr.data.cov <- GenomicRanges::coverage(bw_score, weight = "score")
    seqlevels(region) <- names(gr.data.cov)
    region <- binnedAverage(region, gr.data.cov, "bwscore")
    outerBG <- round(x = mean(region$bwscore), digits = 2)


    # end_time <- Sys.time()
    # message("time: ", end_time - start_time)

    return(outerBG)
}


orderby<-function(peaks,orderFactor){
  if(orderFactor=='distance'){
    message('Ordering by distance to nearest TetO')
    peaks_tmp<-peaks
    if(is.na(peaks_tmp$distance)){
      peaks_tmp$distance<-1e9 #placeholder
    }
    o = order(mcols(peaks_tmp[,'distance']), decreasing = TRUE)
  }else{
    o = order(mcols(peaks[,paste0(orderFactor,'_cov')]), decreasing = TRUE)
  }
  peaks<-peaks[o]
  return(peaks) 
}

filterOverlapPeaks <- function(peaks){
  message('Overlapping peaks identified:', length(unique(findOverlaps(peaks, drop.self=TRUE)@from)))
  #first merge peaks
  peaks_mrgd<-GenomicRanges::reduce(peaks)
  
  if(max(peaks$score)==0){
    message('No score present in peak file, merge peaks')
    peaks<-peaks_mrgd
  }else{
    peaks_hitlist <- as(findOverlaps(query=peaks_mrgd, subject=peaks), "List")
    idx0 <- as(which.max(extractList(peaks$signalValue, peaks_hitlist)), "List")
    idx1 <- unique(unlist(extractList(seq_along(peaks), peaks_hitlist)[idx0]))
    peaks<-peaks[idx1]
    message('Overlapping peaks identified:', length(unique(findOverlaps(peaks, drop.self=TRUE)@from)))  
  }
  
  return(peaks)
}


GetTornadoCov_V2<-function(region, ChIPdf, orderCov=ChIPdf$name, orderby='FLAG_ON', zoom=5000, binWidth=50, Matrix=TRUE, makeList=FALSE, bincoverage=TRUE, prefix){
  
  debug=FALSE
  if(debug){
    
    region=flag_opt_mrgdbyhighestPeak
    ChIPdf=test_ChIPdf
    orderCov=ChIPdf$name
    orderby='C17_FLAG_ON'
    zoom=5000
    binWidth=50
    Matrix=FALSE
    makeList=FALSE
    bincoverage=TRUE
    orderCovname=orderCov[1]
    
    
  }
  
  stopifnot(is(region, "GRanges"))
  stopifnot(is(ChIPdf, "data.frame"))
  
  infoDF<-data.frame(region=deparse(substitute(region)),
                     orderCov=paste(ChIPdf$name,collapse  = '_'),
                     orderby=ifelse(is.null(orderby),'Not ordered',orderby),
                     zoom=zoom,
                     binWidth=binWidth,
                     Matrix=Matrix,
                     makeList=makeList,
                     bincoverage=bincoverage
  )
  
  #orderCov is now just a long string of all the ChIPdf names
  
  if(length(orderCov)>0){
    for(orderCovname in orderCov){
      orderChIPbw<-paste0(prefix,ChIPdf[ChIPdf$name %in% orderCovname,'bigwig'])
      message('Get coverage for ',orderCovname)
      message('bw:',orderChIPbw)
      bw_score<-import(orderChIPbw, selection=BigWigSelection(region))
      gr.data.cov <- GenomicRanges::coverage(bw_score, weight="score")
      seqlevels(region) <- names(gr.data.cov)
      region<-binnedAverage(region, gr.data.cov, paste0(orderCovname,"_cov"))
    }
  }
  
  if(!is.null(orderby)){
    #Check whether coverage is already calculated for order factor otherwise get the coverage
    orderChIPbw<-paste0(prefix,ChIPdf[ChIPdf$name==orderby,'bigwig'])
    message('Order by ',orderby)
    message('bw:',orderChIPbw)
    
    if(!paste0(orderby,"_cov") %in% colnames(mcols(region))){
      message('Get coverage for ',orderby)
      bw_score<-import(orderChIPbw, selection=BigWigSelection(region))
      gr.data.cov <- GenomicRanges::coverage(bw_score, weight="score")
      seqlevels(region) <- names(gr.data.cov)
      region<-binnedAverage(region, gr.data.cov, paste0(orderby,"_cov"))  
    }
    
    o = order(mcols(region[,paste0(orderby,'_cov')]), decreasing = TRUE)
    region<-region[o]
  }
  
  
  #Extend genomic region and tile. This will result in a GRangesList
  message('Tile regions')
  ROI<-GenomicRanges::resize(region,  width = zoom, fix = 'center')
  ROI_tiles<-GenomicRanges::tile(x=ROI, width=binWidth)
  tiles.gr <- unlist(ROI_tiles)
  
  
  if(Matrix){
    Mat_list<-list()
  }
  
  for(a in 1:nrow(ChIPdf)){
    message('Get coverage for ',ChIPdf$name[a],' from ', ChIPdf$bigwig[a])
    bw_score<-import(paste0(prefix,ChIPdf$bigwig[a]), selection=BigWigSelection(tiles.gr))
    
    if(bincoverage){
      message('Calculating averge bin coverage')
      gr.data.cov <- GenomicRanges::coverage(bw_score, weight="score")
      seqlevels(tiles.gr) <- names(gr.data.cov)
      tiles.gr<-binnedAverage(tiles.gr, gr.data.cov, paste0(ChIPdf$name[a],"_cov"))
    }else{
      message('Calculating average bin score')
      hits <- findOverlaps(tiles.gr, bw_score)
      agg <- aggregate(bw_score, hits, score=mean(score))
      mcols(tiles.gr)[paste0(ChIPdf$name[a],"_cov")] <- agg$score
    }
    if(Matrix){
      Mat_list[[ChIPdf$prot[a]]][[ChIPdf$condition[a]]]<-matrix(data = mcols(tiles.gr)[,paste0(ChIPdf$name[a],'_cov')], nrow=length(region), byrow = TRUE)
    }
    
  }
  
  
  if(Matrix){
    return(list(peaks=region, Mat=Mat_list, ChIPdf=ChIPdf, infoDF=infoDF))
  }else{
    
    if(makeList){
      tiles.gr<-relist(tiles.gr, ROI_tiles)
      return(list(peaks=region, tiles=tiles.gr, ChIPdf=ChIPdf, infoDF=infoDF))
    }else{
      return(list(peaks=region, tiles=tiles.gr, ChIPdf=ChIPdf,infoDF=infoDF))
    }
    
  }
  
  
  
  
}

rotate <- function(x) t(apply(x, 2, rev))


avgPlot_V2<-function(Plotregion_cov, 
                  plotOrder=NULL,
                  plotGroups=NULL,
                  orderFactor=NULL,
                  norm='meanBG',
                  normGroup='>3MB',
                  PeaksignalValue=0,
                  plotWidth=2500,
                  FlipRvStrand=FALSE,
                  middleLine=TRUE){
  
  
  debug=FALSE
  if(debug){
    
    Plotregion_cov=toRnadoCov_SMC1_opt_mrgdbyhighestPeak
    plotOrder=plotOrder_1F
    plotGroups=NULL
    orderFactor='C17_SMC1_ON'
    norm='meanBG'
    normGroup='>3MB'
    PeaksignalValue=35
    plotWidth=2500
    FlipRvStrand=FALSE
    middleLine=FALSE
    
  }
  
  #Filter PlotPeaks for signalValue and order based on orderFactor coverage
  PlotPeaks<-Plotregion_cov[['peaks']]
  if(!is.null(orderFactor)){
    message("Order peaks by:",orderFactor)
        
    if(orderFactor=='distance'){
      message('Ordering by distance to nearest TetO')
      PlotPeaks$distance[is.na(PlotPeaks$distance)]<-1e9 #placeholder

      o = order(mcols(PlotPeaks[,'distance']), decreasing = TRUE)
    }else{
      o = order(mcols(PlotPeaks[,paste0(orderFactor,'_cov')]), decreasing = TRUE)
    }
    PlotPeaks<-PlotPeaks[o]
  }
  
  if(PeaksignalValue>0){
    PlotPeaks_sig<-PlotPeaks[PlotPeaks$signalValue>PeaksignalValue]
  }else{
    PlotPeaks_sig<-PlotPeaks
  }
  
  
  ROI_tiles_cov<-Plotregion_cov[['tiles']]
  infoDF<-Plotregion_cov[['infoDF']]
  ChIPdf<-Plotregion_cov[['ChIPdf']]
  
  
  zoom<-unique(infoDF$zoom)/2
  if(!plotWidth){
    plotZoom<-zoom
  }else{
    plotZoom<-plotWidth
  }
  minBin<-(zoom-plotZoom)/infoDF$binWidth + 1
  maxBin<-(zoom+plotZoom)/infoDF$binWidth
  message('plotZoom:',plotZoom)
  
  
  if(!is.null(plotOrder)){
    message('##### reorder ChIPdf and filter based on plotOrder')
    ChIPdf<-ChIPdf[match(plotOrder,ChIPdf$name),]
  }
  

  groupDF<-as.data.frame(table(PlotPeaks_sig$group))
  if(is.null(plotGroups)){
    plotGroups<-groupDF$Var1
  }else{
    plotGroups<-as.factor(plotGroups)
  }
  #select for plotGroups and reorder data.frame
  groupDF <- groupDF[match(plotGroups, groupDF$Var1),]
  
 
  
  
  colMeanDF<-data.frame()
  for(a in 1:nrow(ChIPdf)){
    
    expName<-ChIPdf$name[a]
    ChIPname_prot=ChIPdf[a,'prot']
    
    
    message(expName)
    
    Fullmat_plotOrder<-matrix(data = mcols(ROI_tiles_cov)[,paste0(ChIPdf$name[a],'_cov')], nrow=length(PlotPeaks), byrow = TRUE)

    if(!is.null(orderFactor)){
      message('Order matrix')
      Fullmat_plotOrder<-Fullmat_plotOrder[o,]
    }
    
    if(PeaksignalValue>0){
      Fullmat_plotOrder<-Fullmat_plotOrder[which(PlotPeaks$signalValue>PeaksignalValue),]  
    }else{
      Fullmat_plotOrder<-Fullmat_plotOrder
    }
    
    
    if(FlipRvStrand){
      message('downstream motifs will be flipped > upstream')
      Fullmat_plotOrder[which(strand(PlotPeaks_sig)=='-'),]<-t(apply(Fullmat_plotOrder[which(strand(PlotPeaks_sig)=='-'),], 1, rev))
    }
    
    if(norm=='meanBG'){
      message('Normalize based on meanBG')
      Fullmat_plotOrder<-Fullmat_plotOrder/ChIPdf$meanBG[a]
    }
    
    if(norm=='outerBG'){
      message('Normalize based on outerBG')
      Fullmat_plotOrder<-Fullmat_plotOrder/ChIPdf$outerBG[a]
    }
    
    if(norm=='outerBins'){
      message('Normalize based on outer 100 bins selected rows')
      groupIDX<-which(PlotPeaks_sig$group %in% normGroup)
      maxOuterBin=ncol(Fullmat_plotOrder)
      minOuterBin=ncol(Fullmat_plotOrder)-100
      
      meanBG<-mean(Fullmat_plotOrder[groupIDX,c(1:100,minOuterBin:maxOuterBin)])
      
      message('outerBins meanBG:',meanBG)
      Fullmat_plotOrder<-Fullmat_plotOrder/meanBG
      
    }
    
    
    for(plotGroup in plotGroups){

      groupIDX<-which(PlotPeaks_sig$group %in% plotGroup)

      
      colM<-colMeans(Fullmat_plotOrder[groupIDX,])
      colM<-as.data.frame(matrix(data = colM, nrow = 1, ncol = length(colM)))
      
      #subset based on plotWidth
      colM<-colM[,minBin:maxBin]

      newRow<-data.frame(expName=expName, 
                         group=plotGroup, 
                         colM
      )
      col
      
      
      colMeanDF<-rbind(colMeanDF, newRow
      )
      
    }
    
  }
  
  return(colMeanDF)
}


getCov<-function(bw, region, name){
  bw_score<-import(bw, selection=BigWigSelection(region))
  GRcov <- GenomicRanges::coverage(bw_score, weight="score")
  seqlevels(region) <- names(GRcov)
  region<-binnedAverage(region, GRcov, name)
  return(region)
}

makePlotToRnado_V8 <- function(Plotregion_cov,
                               plotOrder = NULL,
                               plotGroups = NULL,
                               orderFactor = NULL,
                               PlotMaxasCutoff = FALSE,
                               norm = 'meanBG',
                               normGroup = NULL,
                               SetCutOff = NULL,
                               cutoffRatio = 1,
                               cutoffGroups = NULL,
                               cutoffLine = TRUE,
                               subSampleGW = TRUE,
                               subSampleGroups = '>3MB',
                               subsampleN = NULL,
                               PeaksignalValue = 0,
                               equalPlotMax = FALSE,
                               setPlotMax = NULL,
                               plotWidth = FALSE,
                               FlipRvStrand = FALSE,
                               BigMatrix = FALSE,
                               middleLine = TRUE,
                               protColors = protCols,
                               TetO_enrichment = TRUE,
                               debugLines = TRUE,
                               PNG = FALSE,
                               PDF = FALSE,
                               useRaster = TRUE,
                               plotName = '',
                               ColWidth = 4,
                               prefix = '',
                               outF = '/TACL/PK/toRnado/figures/',
                               verbose = TRUE) {
  
  #normGroup only required if norm='outerBins'
  
  
  debug <- FALSE
  if (debug) {
    Plotregion_cov = toRnadoCov_MGA_K562
    plotOrder = plotOrder = ChIPdf$name
    plotGroups = c('MGA+;CTCF-', 'MGA+;CTCF+')
    orderFactor = 'MGA'
    PlotMaxasCutoff = TRUE
    norm = 'meanBG'
    normGroup = c('MGA+;CTCF-', 'MGA+;CTCF+')
    SetCutOff = NULL
    cutoffRatio = 0.6
    cutoffGroups = NULL
    cutoffLine = TRUE
    subSampleGW = FALSE
    subSampleGroups = NULL
    subsampleN = NULL
    PeaksignalValue = 0
    equalPlotMax = FALSE
    setPlotMax = NULL
    plotWidth = FALSE
    FlipRvStrand = FALSE
    BigMatrix = TRUE
    middleLine = TRUE
    protColors = NULL
    PNG = FALSE
    PDF = FALSE
    plotName = ''
    ColWidth = 4
    outF = '/home/p.krijger_cbs-niob.local/projects/Sjoerd/toRnado/'
  }
  
  if (!dir.exists(outF)) {
    dir.create(outF)
  }
  
  if (is.null(norm)) {
    message("No normalization")
    norm <- "noNorm"
  }
  
   linecols <- c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000",
    "grey",
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D"
  )
  
  
  if (PNG | PDF) {
    ifelse(!dir.exists(file.path(prefix, outF)),
           dir.create(file.path(prefix, outF)),
           'out folder already exists')
  }
  
  
  #Filter PlotPeaks for signalValue and order based on orderFactor coverage
  PlotPeaks <- Plotregion_cov[['peaks']]
  if (!is.null(orderFactor)) {
    message("Order peaks by:", orderFactor)
    
    
    if (orderFactor == 'distance') {
      message('Ordering by distance to nearest TetO')

      #PlotPeaks <- toRnadoCov$peaks

      if (any(is.na(PlotPeaks$distance))) {
      PlotPeaks[is.na(PlotPeaks$distance)]$distance <- 1e9 #placeholder

      }

      o = order(mcols(PlotPeaks[, 'distance']), decreasing = FALSE)
    } else{
      o = order(mcols(PlotPeaks[, paste0(orderFactor, '_cov')]), decreasing = TRUE)
    }
    PlotPeaks <- PlotPeaks[o]
  }
  
  if (PeaksignalValue > 0) {
    PlotPeaks_sig <- PlotPeaks[PlotPeaks$signalValue > PeaksignalValue]
  } else{
    PlotPeaks_sig <- PlotPeaks
  }
  
  ROI_tiles_cov <- Plotregion_cov[['tiles']]
  infoDF <- Plotregion_cov[['infoDF']]
  ChIPdf <- Plotregion_cov[['ChIPdf']]
  
  
  zoom <- unique(infoDF$zoom) / 2
  if (!plotWidth) {
    plotZoom <- zoom
  } else{
    plotZoom <- plotWidth
  }
  minBin <- (zoom - plotZoom) / infoDF$binWidth + 1
  maxBin <- (zoom + plotZoom) / infoDF$binWidth
  
  message('plotZoom:', plotZoom)
  message('zoom:', zoom)
  message('minBin:', minBin)
  message('maxBin:', maxBin)

  
  if (!is.null(plotOrder)) {
    message('##### reorder ChIPdf and filter based on plotOrder')
    ChIPdf <- ChIPdf[match(plotOrder, ChIPdf$name), ]
  }
  
  
  #the groups are the rows that will be plotted
  message('##### get the groups that will be plotted')
  groupDF <- as.data.frame(table(PlotPeaks_sig$group))
  if (is.null(plotGroups)) {
    plotGroups <- groupDF$Var1
  } else{
    plotGroups <- as.factor(plotGroups)
  }
  #select for plotGroups and reorder data.frame
  groupDF <- groupDF[match(plotGroups, groupDF$Var1), ]
  
  if (is.null(cutoffGroups)) {
    cutoffGroups <- plotGroups
  }
  
  #message(groupDF)
  
  
  message('Calculating plotMax for each group')
  plotMax_df <- data.frame()
  for (a in 1:nrow(ChIPdf)) {
    if (verbose) {
      message(ChIPdf$name[a])
    }
    
    Fullmat_plotMax <-
      matrix(
        data = mcols(ROI_tiles_cov)[, paste0(ChIPdf$name[a], '_cov')],
        nrow = length(PlotPeaks),
        byrow = TRUE
      )
    
    if (!is.null(orderFactor)) {
      message('Order matrix')
      Fullmat_plotMax <- Fullmat_plotMax[o, ]
    }
    
    if (PeaksignalValue > 0) {
      Fullmat_plotMax <-
        Fullmat_plotMax[which(PlotPeaks$signalValue > PeaksignalValue), ]
    }
    
    
    if (FlipRvStrand) {
      message('downstream motifs will be flipped > upstream')
      Fullmat_plotMax[which(strand(PlotPeaks_sig) == '-'), ] <-
        t(apply(Fullmat_plotMax[which(strand(PlotPeaks_sig) == '-'), ], 1, rev))
    }
    
    if (norm == 'meanBG') {
      message('Normalize based on meanBG')
      Fullmat_plotMax <- Fullmat_plotMax / ChIPdf$meanBG[a]
    }
    
    if (norm == 'outerBG') {
      message('Normalize based on outerBG')
      Fullmat_plotMax <- Fullmat_plotMax / ChIPdf$outerBG[a]
    }
    
    
    if (norm == 'outerBins') {
      if (ChIPdf$prot[a] == '4C') {
        message('4C data already normalized')
      } else{
        message('Normalize based on outer 100 bins selected rows')
        groupIDX <- which(PlotPeaks_sig$group %in% normGroup)
        maxOuterBin = ncol(Fullmat_plotMax)
        minOuterBin = ncol(Fullmat_plotMax) - 100
        meanBG <-
          mean(Fullmat_plotMax[groupIDX, c(1:100, minOuterBin:maxOuterBin)])
        message('MeanBG:', meanBG)
        Fullmat_plotMax <- Fullmat_plotMax / meanBG
      }
    }
    
    newRow <- data.frame(
      name = ChIPdf$name[a],
      protein = ChIPdf$prot[a],
      condition = ChIPdf$condition[a],
      Fullmat_max = max(colMeans(Fullmat_plotMax))
    )
    
    
    for (cutoffGroup in cutoffGroups) {
      groupIDX <- which(PlotPeaks_sig$group %in% cutoffGroup)
      newRow[, paste0(cutoffGroup)] <-
        max(colMeans(Fullmat_plotMax[groupIDX, ]))
    }
    plotMax_df <- rbind(plotMax_df, newRow)
    
  }

  aggDF <-
    aggregate(plotMax_df[, c(4:ncol(plotMax_df))], by = list(plotMax_df$protein), max)
  
  
  
  
  #Make Layout based on plotGroups and plotOrder
  
  mat.col <- nrow(ChIPdf) + 1
  layout.name.matrix <-
    matrix(
      rep(1, each = mat.col),
      nrow = 1,
      ncol = mat.col,
      byrow = TRUE
    )
  
  
  
  
  if (!PNG & !PDF) {
    #dev.new(width=(4*mat.col)+2, height=30, unit="cm")
    openWindow <- readline(prompt = "Open new window (y/n)")
    
    if (openWindow %in% c('yes', 'YES', 'Yes', 'Y', 'y')) {
      dev.new(width = 25, height = 11.5)
    }
    
  }
  
  if (PNG) {
    message('Making PNG')
    peakname <- infoDF$region
    groupNames <- paste(plotGroups, collapse = '_')
    groupNames <-
      gsub(pattern = '<',
           replacement = 'smallerthan',
           x = groupNames)
    groupNames <-
      gsub(pattern = '>',
           replacement = 'largerthan',
           x = groupNames)
    groupNames <-
      gsub(pattern = ' ',
           replacement = '',
           x = groupNames)
    
    normName <- ifelse(is.null(norm), 'noNorm', norm)
    outfileName <-
      paste0(
        plotName,
        '_',
        peakname,
        '_by_',
        orderFactor,
        '_',
        cutoffRatio,
        '_Flip',
        FlipRvStrand,
        '_',
        normName,
        '.png'
      )
    
    
    outfile <- paste0(prefix, outF, outfileName)
    
    if (file.exists(outfile)) {
      message('Warning: file already exists. modifying file name')
      
      outfile_new <-
        sub(
          pattern = '.png',
          replacement = paste0(Sys.Date(), '.png'),
          x = outfile
        )
      
      message('file name:', outfile_new)
      
      if (file.exists(outfile_new)) {
        datetime_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
        outfile_new <-
          sub(
            pattern = '.png',
            replacement = paste0(datetime_str, '.png'),
            x = outfile
          )
        message('new name:', outfile_new)
        outfile <- outfile_new
        
      } else{
        message('new name:', outfile_new)
        outfile <- outfile_new
      }
      
    } else{
      message('pdf file:', outfile)
    }
    
    png(
      outfile,
      width =  (mat.col * 4) / (2.54 / 96) + 200,
      height = 22 / (2.54 / 96) + 200,
      res = 100
    )
    
    
    
  }
  
  
  if (PDF) {
    message('Making PDF')
    peakname <- infoDF$region
    groupNames <- paste(plotGroups, collapse = '_')
    groupNames <-
      gsub(pattern = '<',
           replacement = 'smallerthan',
           x = groupNames)
    groupNames <-
      gsub(pattern = '>',
           replacement = 'largerthan',
           x = groupNames)
    groupNames <-
      gsub(pattern = ' ',
           replacement = '',
           x = groupNames)
    
    normName <- ifelse(is.null(norm), 'noNorm', norm)
    outfileName <-
      paste0(
        plotName,
        '_',
        peakname,
        '_by_',
        orderFactor,
        '_',
        cutoffRatio,
        '_Flip',
        FlipRvStrand,
        '_',
        normName,
        '.pdf'
      )
    
    
    outfile <- paste0(prefix, outF, outfileName)
    
    if (file.exists(outfile)) {
      message('Warning: file already exists. modifying file name')
      
      outfile_new <-
        sub(
          pattern = '.pdf',
          replacement = paste0(Sys.Date(), '.pdf'),
          x = outfile
        )
      
      message('file name:', outfile_new)
      
      if (file.exists(outfile_new)) {
        datetime_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
        outfile_new <-
          sub(
            pattern = '.pdf',
            replacement = paste0(datetime_str, '.pdf'),
            x = outfile
          )
        message('new name:', outfile_new)
        outfile <- outfile_new
        
      } else{
        message('new name:', outfile_new)
        outfile <- outfile_new
      }
      
    }
    
    pdf(
      outfile,
      width =  (mat.col * 4) / (2.54) + 2.54,
      height = 23 / (2.54) + 2.54
    )
  }
  
  
  
  
  if (BigMatrix) {
    if (verbose) {
      message('Plot 1 big matrix')
    }
    
    mat.row <- 5
    layout.plot.matrix <-
      matrix(
        2:(4 * mat.col + 1),
        nrow = mat.row - 1,
        ncol = mat.col,
        byrow = FALSE
      )
    layout.matrix <- rbind(layout.name.matrix, layout.plot.matrix)
    layout(
      mat = layout.matrix
      ,
      height = c(lcm(1.5), lcm(4.5), lcm(15), lcm(1), lcm(1))
      ,
      width = c(lcm(1), rep(lcm(ColWidth), mat.col - 1))
    )
    
  } else{
    message('Split up Matrix in seperate blocks: NOT FINISHED')
    #use groupDF
    
    # mat.row<-length(plotGroups)+2
    # layout.plot.matrix <- matrix(2:(mat.row*mat.col), nrow = mat.row, ncol = mat.col, byrow = FALSE)
    # layout.matrix <-rbind(layout.name.matrix,layout.plot.matrix)
    # layout(mat = layout.matrix
    #        ,height = c(lcm(1),lcm(4),rep(lcm(5),length(plotGroups)), lcm(0.5))
    #        ,width = rep(lcm(ColWidth),mat.col)
    # )
  }
  
  
  
  
  #Title
  
  if (verbose) {
    message('Plot Title')
  }
  
  par(mar = c(0, 0, 0, 0))
  plot(
    c(0, 1),
    c(0, 3),
    ann = F,
    bty = 'n',
    type = 'n',
    xaxt = 'n',
    yaxt = 'n'
  )
  
  text(
    x = 0.5,
    y = 2.5,
    paste(
      "Peak:",
      infoDF$region,
      "order:",
      orderFactor,
      "Peak signalValue:",
      PeaksignalValue,
      "cutoffRatio:",
      cutoffRatio,
      "norm:",
      norm
    ),
    cex = 1.4,
    col = "black"
  )
  
  #text(x = 0.5, y = 0.5, paste(plotGroups,groupDF[plotGroups,'Freq'], sep= ':', collapse = ' '), cex = 1, col = "black")
  
  # for(a in length(plotGroups):1){
  #
  #   text(x = 0.5, y = 0.5, paste(plotGroups[1:a],groupDF[plotGroups[1:a],'Freq'], sep= ':', collapse = ' '), cex = 1, col = "black")
  # }
  
  
  
  
  # legend("center"
  #        ,legend=paste(groupDF$Var1,groupDF$Freq, sep= ':')
  #        ,fill = linecols[1:length(plotGroups)]
  #        ,bty ="n", pch=NA, horiz=TRUE)
  
  
  
  legend_text <- paste(groupDF$Var1, groupDF$Freq, sep = ':')
  
  if (length(legend_text) > 4) {
    n_col <- ceiling(length(legend_text) / 2)
    MyOrder = matrix(
      1:length(legend_text),
      nrow = 2,
      ncol = n_col,
      byrow = T
    )
    
    legend(
      "bottom",
      legend_text[MyOrder],
      fill = linecols[MyOrder],
      ncol = n_col,
      border = NA ,
      bty = "n",
      pch = NA
    )
    
    # legend("center"
    #        ,legend=legend_text
    #        ,fill = linecols[1:length(plotGroups)]
    #        ,bty ="n", pch=NA, ncol=n_col)
    
    
    
    
  } else{
    legend(
      "center"
      ,
      legend = legend_text
      ,
      fill = linecols[1:length(plotGroups)]
      ,
      bty = "n",
      pch = NA,
      horiz = TRUE
    )
  }
  
  
  
  
  
  
  
  #first plot the group legend
  
  
  
  plot.new() # empty row between title and heatmaps. for the factors this is where the avg plot is.
  
  
  
  
  
  if (BigMatrix) {
    if (verbose) {
      message('Plot Big Matrix')
    }
    
    
    par(mar = c(0, 2, 0, 0) + 0.1)  #c(bottom, left, top, right).
    
    
    groupList <- list()
    for (plotGroup in plotGroups) {
      groupList[[plotGroup]] <- which(PlotPeaks_sig$group %in% plotGroup)
    }
    
    
    if (subSampleGW) {
      if (!any(subSampleGroups %in% groupDF$Var1)) {
        message('SampleGroups not (all) in plotGroups')
        
        
        message(subSampleGroups %in% groupDF$Var1)
        message(subSampleGroups)
        message(groupDF$Var1)
        
        
      } else{
        message('subsample')
        
        if (is.null(subsampleN)) {
          if (nrow(groupDF) == 1) {
            message("WARNING: only 1 group defined. Subsetting by lowest group not possible")
          } else{
            nMax <- max(groupDF[!groupDF$Var1 %in% subSampleGroups, ]$Freq)
            message('nMax:', nMax)
            
            for (subSampleGroup in subSampleGroups) {
              message('subSampleGroup:', subSampleGroup)
              subsetID <-
                floor(seq(
                  from = 1,
                  to = groupDF[groupDF$Var1 %in% subSampleGroup, ]$Freq,
                  length.out = nMax
                ))
              groupList[[subSampleGroup]] <-
                groupList[[subSampleGroup]][subsetID]
              
            }
            
            
          }
          
          
        } else{
          nMax <- subsampleN
          
          for (subSampleGroup in subSampleGroups) {
            message('subSampleGroup:', subSampleGroup)
            subsetID <-
              floor(seq(
                from = 1,
                to = groupDF[groupDF$Var1 %in% subSampleGroup, ]$Freq,
                length.out = nMax
              ))
            groupList[[subSampleGroup]] <-
              groupList[[subSampleGroup]][subsetID]
            
          }
        }
        
      }
      
    }
    
    
    
    message('Making legend Matrix')
    
    #we here add an empty (white) row
    
    groupCat = 0
    for (plotGroup in plotGroups) {
      message(plotGroup)
      groupCat = groupCat + 1
      if (groupCat == 1) {
        message(length(groupList[[plotGroup]]))
        #matTMP<-matrix(data = rep(x = linecols[groupCat], length(groupList[[plotGroup]])),nrow = length(groupList[[plotGroup]]), ncol = 1)
        matTMP <-
          matrix(
            data = rep(x = groupCat, length(groupList[[plotGroup]])),
            nrow = length(groupList[[plotGroup]]),
            ncol = 10
          )
      } else{
        message('Add empty row and new rows')
        message(length(groupList[[plotGroup]]))
        matTMP <- rbind(
          matTMP,
          #matrix(data = '#ffffff', nrow = 1, byrow = TRUE),
          #matrix(data = rep(x = linecols[groupCat], length(groupList[[plotGroup]])),nrow = length(groupList[[plotGroup]]), ncol = 1)
          matrix(
            data = 17, #white
            nrow = 1,
            ncol = 10
          ),
          matrix(
            data = rep(x = groupCat, length(groupList[[plotGroup]])),
            nrow = length(groupList[[plotGroup]]),
            ncol = 10
          )
          
        )
      }
      
      
      
    }
    
    
    
    # linecols2 <-
    #   c(
    #     "#E69F00",
    #     "#56B4E9",
    #     "#009E73",
    #     "#F0E442",
    #     "#0072B2",
    #     "#D55E00",
    #     "#CC79A7",
    #     "#000000",
    #     "grey",
    #     "grey",
    #     '#ffffff'
    #   )
    #plot(x = c(1:11), y=c(1:11), col=linecols2)
    
       linecols2 <- c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000",
    "grey",
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D",
    '#ffffff'
  )


    matTMP <- rotate(matTMP)
    
    #dim(matTMP)
    
    #image(matTMP, col = linecols2, axes=FALSE)
    
    if (useRaster) {
      image(matTMP,
            col = linecols2,
            axes = FALSE,
            useRaster = TRUE)
    } else{
      image(matTMP,
            col = linecols2,
            axes = FALSE,
            useRaster = FALSE)
    }
    
    
    
    
    
    
    #matTMP2<-matTMP[1:10,c(1:100,33485:33584)]
    #image(matTMP2,col = linecols2, axes=FALSE, useRaster=TRUE)
    #image(1:nrow(matTMP2), 1:ncol(matTMP2), matTMP2, col = linecols2, axes = FALSE)
    
    
    #my_matrix <- matrix(1:25, nrow = 5)
    #image(my_matrix, col=linecols2)
    #image(my_matrix, col=linecols2, useRaster = TRUE)
    
    
    #image(1:nrow(matTMP2), 1:ncol(matTMP2), matTMP2, col = terrain.colors(60), axes = FALSE)
    
    # matTMP2 <- matrix(c(rep(1, 50),rep(2, 50)), nrow = 10)  # Example matrix data
    # image(t(matTMP2), col = c("white", "black","red"), axes=FALSE)
    #
    #
    # image(matTMP2, col = linecols2, axes=FALSE)
    
    
    # m <- matrix(1:30, ncol=6)
    # colnames(m) <- paste("C", 1:6, sep="")
    # rownames(m) <- paste("R", 1:5, sep="")
    # m
    #
    # image(1:ncol(m), 1:nrow(m), t(m), col = terrain.colors(60), axes = FALSE)
    # axis(1, 1:ncol(m), colnames(m))
    # axis(2, 1:nrow(m), rownames(m))
    # for (x in 1:ncol(m))
    #   for (y in 1:nrow(m))
    #     text(x, y, m[y,x])
    
    
    
  } else{
    #no Big Matrix
    plot.new()
  }
  
  #why do i do this? I guess to have a black space in between the legend matrix and the data?
  plot.new()
  plot.new()
  
  
  
  
  
  
  
  
  
  #Read matrix per Chipfactor
  for (a in 1:nrow(ChIPdf)) {
    ChIPname_prot = ChIPdf[a, 'prot']
    aggDF_prot <- aggDF[aggDF$Group.1 == ChIPname_prot, ]
    Fullmat_plotOrder <-
      matrix(
        data = mcols(ROI_tiles_cov)[, paste0(ChIPdf$name[a], '_cov')],
        nrow = length(PlotPeaks),
        byrow = TRUE
      )
    if (!is.null(orderFactor)) {
      message('Order matrix')
      Fullmat_plotOrder <- Fullmat_plotOrder[o, ]
    }
    
    if (PeaksignalValue > 0) {
      Fullmat_plotOrder <-
        Fullmat_plotOrder[which(PlotPeaks$signalValue > PeaksignalValue), ]
    }
    
    
    if (FlipRvStrand) {
      message('downstream motifs will be flipped > upstream')
      Fullmat_plotOrder[which(strand(PlotPeaks_sig) == '-'), ] <-
        t(apply(Fullmat_plotOrder[which(strand(PlotPeaks_sig) == '-'), ], 1, rev))
    }
    
    if (norm == 'meanBG') {
      message('Normalize based on meanBG')
      Fullmat_plotOrder <- Fullmat_plotOrder / ChIPdf$meanBG[a]
    }
    
    if (norm == 'outerBG') {
      message('Normalize based on outerBG')
      Fullmat_plotOrder <- Fullmat_plotOrder / ChIPdf$outerBG[a]
    }
    
    
    if (norm == 'outerBins') {
      if (ChIPdf$prot[a] == '4C') {
        message('4C data already normalized')
      } else{
        message('Normalize based on outer 100 bins selected rows')
        groupIDX <- which(PlotPeaks_sig$group %in% normGroup)
        maxOuterBin = ncol(Fullmat_plotOrder)
        minOuterBin = ncol(Fullmat_plotOrder) - 100
        meanBG <-
          mean(Fullmat_plotOrder[groupIDX, c(1:100, minOuterBin:maxOuterBin)])
        message('MeanBG:', meanBG)
        Fullmat_plotOrder <- Fullmat_plotOrder / meanBG
      }
    }
    
    #Get plotMax
    if (!is.null(setPlotMax)) {
      message('plotMax set to: ', setPlotMax)
      plotMax <- setPlotMax
    } else{
      if (equalPlotMax) {
        message('Select PlotMax as max for all groups and all proteins')
        normIDX <- which(colnames(aggDF) %in% cutoffGroups)
        #plotMax<-ceiling(max(aggDF[,normIDX])*1.1)
        plotMax <- round(max(aggDF[, normIDX]) * 1.1, 1)
        
        
        
      } else{
        message('Select PlotMax as max for all groups per protein')
        normIDX <- which(colnames(aggDF_prot) %in% cutoffGroups)
        #plotMax<-ceiling(max(aggDF_prot[,normIDX])*1.1)
        plotMax <- round(max(aggDF_prot[, normIDX]) * 1.1, 1)
      }
      
    }
    
    
    if (is.null(SetCutOff)) {
      if (PlotMaxasCutoff) {
        message('Using PlotMaxasCutoff')
        normIDX <- which(colnames(aggDF) %in% cutoffGroups)
        cutoff <- max(aggDF_prot[, normIDX]) * cutoffRatio
      } else{
        message('Using aggDF_prot$Fullmat_max*cutoffRatio')
        cutoff <- aggDF_prot$Fullmat_max * cutoffRatio
      }
    } else{
      if (!cutoffRatio == 1) {
        message('WARNING: CUT off ration used in combination with cutoff')
      }
      
      
      cutoff <-
        SetCutOff[a] * cutoffRatio  #Do we really want to use a cutoffRatio here??
    }
    message('Cutoff:', cutoff)
    
    
    
    #Mean plot
    message('Make density plot with plotMax: ', plotMax)
    #par(mar=c(0, 2, 2, 0) +0.1)  #c(bottom, left, top, right).
    
    par(mar = c(0, 2, 5, 0) + 0.1)  #c(bottom, left, top, right).
    
    plotTitle <- ChIPdf$name[a]
    
    
    # name      RH_cell cell condition    prot
    
    # if(plotTitle=='C17_TmCherry_MAU2-AID-GFP_GFP'){
    #   plotTitle<-'MAU2-GFP'
    # }
    # if(plotTitle=='C17_TmCherry_V5-MAU2_V5'){
    #   plotTitle<-'V5-MAU2_OFF'
    # }
    # if(plotTitle=='C17_TmCherry_V5-MAU2_NIPBL'){
    #   plotTitle<-'V5-MAU2_OFF'
    # }
    # if(plotTitle=='C17_TMAU2_V5-MAU2_V5'){
    #
    # }
    # if(plotTitle=='C17_TMAU2_V5-MAU2_NIPBL'){
    #
    # }
    #
    
    plotTitle <-
      gsub(pattern = 'C17_',
           replacement = '',
           x = plotTitle)
    plotTitle <-
      gsub(pattern = 'RAD21AID_SC1_',
           replacement = '',
           x = plotTitle)
    plotTitle <-
      gsub(pattern = 'RAD21AID_SC6_',
           replacement = '',
           x = plotTitle)
    
    #plot(0, type='n',  xlim=c(minBin,maxBin), ylim=c(0,plotMax), axes = FALSE, xlab = "", ylab = "mean cov", main = plotTitle)
    
    
    
    plotcell <- ChIPdf$cell[a]
    plotcond <- ChIPdf$condition[a]
    plotprot <- ChIPdf$prot[a]
    plotrep <- ChIPdf$rep[a]


    message('minBin:', minBin)
    message('maxBin:', maxBin)
    
    plot(
      0,
      type = 'n',
      xlim = c(minBin, maxBin),
      ylim = c(0, plotMax),
      axes = FALSE,
      xlab = "",
      ylab = "mean cov",
      main = "",
      cex.main = 1.2,
      font.main = 2
    )
    mtext(
      paste(plotprot,plotrep),
      side = 3,
      line = 1.5,
      cex = 0.8,
      font = 2
    )

    mtext(
      plotcond,
      side = 3,
      line = 2.5,
      cex = 0.8,
      font = 2
    )
    mtext(
      plotcell,
      side = 3,
      line = 3.5,
      cex = 0.8,
      font = 2
    )

    if(TetO_enrichment){
        tetO <- round(ChIPdf$TetO_enrichment[a],1)

    mtext(
      tetO,
      side = 3,
      line = 4.5,
      cex = 0.8,
      font = 2,
      col = ifelse(is.na(tetO), 'black', ifelse(tetO >= 10, 'darkgreen', 'darkred'))
    )

    }

    


    
    #axis(2, at = c(round(cutoff,2),seq(0,plotMax, by = 1)))
    
    if(debugLines){
      message('Line 1444')
    }
    

    if (cutoffLine) {
      axis(2, at = c(0, round(cutoff, 1), plotMax))
    } else{
      axis(2, at = c(0, plotMax))
    }
    
    
    #axis(2, at = c(round(cutoff,1),plotMax))
    
    
    #https://stats.stackexchange.com/questions/118033/best-series-of-colors-to-use-for-differentiating-series-in-publication-quality
    
    for (b in 1:length(plotGroups)) {
      groupIDX <- which(PlotPeaks_sig$group %in% plotGroups[b])
      lines(colMeans(Fullmat_plotOrder[groupIDX, ]),
            col = linecols[b],
            type = 'l')
    }
    
    if (cutoffLine) {
      abline(h = cutoff, lty = 3)
    }
    
    
    # linecols<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
    #             "#0072B2", "#D55E00", "#CC79A7", "#000000")
    # legend("topleft"
    #        ,legend=plotGroups
    #        ,fill = linecols[1:length(plotGroups)]
    #        ,bty ="n", pch=NA)
    
    if (middleLine) {
      abline(v = (maxBin + minBin) / 2, lty = 2)
    }
    
    
    
    
    
    #heatmap
    message('Make tornado')
    # protColors <- data.frame(prot=c('FLAG','SMC1','NIPBL','CTCF','V5','H3K27ac','H3K4me3','RAD21','MAU2','H3K27me3','PolII','STAG1', 'STAG2'),
    #                          color=c('#377eb8','#e41a1c','#4daf4a','#808080','#ff7f00','#8B8000','#a65628','#AA336A','#3D426B','#72d6c9','#405983','#636363','#636363'),
    #                          bgColour=c('#b3cde3','#fbb4ae','#ccebc5','#F2F3F4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#b3cde3','#ffdd83','#e3f8ff','#f0f0f0','#f0f0f0'))
    #
    # bgColour<-protColors[protColors$prot==ChIPname_prot,'bgColour']
    # protColour<-protColors[protColors$prot==ChIPname_prot,'color']
    
    if (is.null(protColors)) {
      bgColour <- '#f0f0f0'
      protColour <- '#636363'
    } else{
      message('Get prot Colour from table for protein: ', ChIPname_prot)
      bgColour <-
        protColors[protColors$protein == ChIPname_prot, 'bgColour']
      protColour <-
        protColors[protColors$protein == ChIPname_prot, 'color']
    }
    
    #protCols[protCols$protein == ChIPname_prot, 'bgColour'] #commented out 31.05.24
    
    if(plotcond=="ON-OFF"){
        #delta Cov can be negative
        paletteLength <- 101
        myColor <- c(colorRampPalette(c("red", "white"))(50), colorRampPalette(c("white", "blue"))(50))
        myBreaks <- seq(-cutoff, cutoff, length.out = paletteLength)
        hic.mat.z <- Fullmat_plotOrder
        hic.mat.z[hic.mat.z > cutoff] <- cutoff
        minCutoff <- -cutoff
        hic.mat.z[hic.mat.z < minCutoff] <- minCutoff

        #heatmap(delta_mat, scale = "none", Rowv = NA, Colv = NA, col = colorRampPalette(c("blue", "white", "red"))(100), margins = c(5, 10), labRow = "", labCol = "", main = "dMAU2")

    }else{
        paletteLength <- 100 #20
        myColor <- colorRampPalette(c(bgColour, protColour))(paletteLength)
        myBreaks <- seq(0, cutoff, length.out = paletteLength)
        
        hic.mat.z <- Fullmat_plotOrder
        hic.mat.z[hic.mat.z > cutoff] <- cutoff
        minCutoff <- 0



    }



    

    

    


    

    
    #which row is na
    #table(is.na(hic.mat.z)) #none = NA
    #hic.mat.z[rowSums(is.na(hic.mat.z)) != ncol(hic.mat.z), ]
    
    #na.omit is that it removes rows with any NA
    #rowSums(is.na(hic.mat.z)) = ncol(hic.mat.z)
    
    
    #ind <- apply(hic.mat.z, 1, function(x) all(is.na(x)))
    #X <- X[ !ind, ]
    
    
    if (BigMatrix) {
      message('Plot Big Matrix')
      par(mar = c(0, 2, 0, 0) + 0.1)  #c(bottom, left, top, right).
      
      
      groupList <- list()
      for (plotGroup in plotGroups) {
        groupList[[plotGroup]] <- which(PlotPeaks_sig$group %in% plotGroup)
      }
      
      
      if (subSampleGW) {
        if (!any(subSampleGroups %in% groupDF$Var1)) {
          message('SampleGroups not (all) in plotGroups')
          
          
          message(subSampleGroups %in% groupDF$Var1)
          message(subSampleGroups)
          message(groupDF$Var1)
          
          
        } else{
          message('subsample')
          
          if (is.null(subsampleN)) {
            if (nrow(groupDF) == 1) {
              message("WARNING: only 1 group defined. Subsetting by lowest group not possible")
            } else{
              nMax <- max(groupDF[!groupDF$Var1 %in% subSampleGroups, ]$Freq)
              message('nMax:', nMax)
              
              for (subSampleGroup in subSampleGroups) {
                message('subSampleGroup:', subSampleGroup)
                subsetID <-
                  floor(seq(
                    from = 1,
                    to = groupDF[groupDF$Var1 %in% subSampleGroup, ]$Freq,
                    length.out = nMax
                  ))
                groupList[[subSampleGroup]] <-
                  groupList[[subSampleGroup]][subsetID]
                
              }
              
              
            }
            
            
          } else{
            nMax <- subsampleN
            
            for (subSampleGroup in subSampleGroups) {
              message('subSampleGroup:', subSampleGroup)
              subsetID <-
                floor(seq(
                  from = 1,
                  to = groupDF[groupDF$Var1 %in% subSampleGroup, ]$Freq,
                  length.out = nMax
                ))
              groupList[[subSampleGroup]] <-
                groupList[[subSampleGroup]][subsetID]
              
            }
          }
          
        }
        
      }
      
      
      
      
      message('#combine matrixes')
      emptyrow <-
        matrix(data = rep(0, ncol(hic.mat.z)),
               nrow = 1,
               byrow = TRUE)
      
      
      
      
      groupCat = 0
      for (plotGroup in plotGroups) {
        message(plotGroup)
        groupCat = groupCat + 1
        if (groupCat == 1) {
          matTMP <- hic.mat.z[groupList[[plotGroup]], minBin:maxBin]
        } else{
          message('Add empty row and new rows')
          matTMP <- rbind(matTMP,
                          emptyrow[, minBin:maxBin],
                          hic.mat.z[groupList[[plotGroup]], minBin:maxBin])
        }
        
        
        
      }
      
      
      
      hic.mat.z <- rotate(matTMP)
      #juicerCol <-c("white","red")
      #wr <- colorRampPalette(juicerCol)
      #image(hic.mat.z, col=wr(1e4), axes=FALSE, zlim=c(0,cutoff))
      
      #define the rows?
      #image(bin_centers, 1:nrow(scaled_data), t(scaled_data), xlab = "Distance to TetO (Mb)", ylab = "TACL TetOs", main = "#FLAG peaks, surrounding TACL TetOs", col = aminCol, xaxt = "n")
      
      if (useRaster) {
        image(
          hic.mat.z,
          col = myColor,
          axes = FALSE,
          zlim = c(minCutoff, cutoff),
          useRaster = TRUE
        )
      } else{
        image(
          hic.mat.z,
          col = myColor,
          axes = FALSE,
          zlim = c(minCutoff, cutoff),
          useRaster = FALSE
        )
      }
      
      box(lwd = 2)
      
      #draw line
      y = 0
      
      if (length(plotGroups) > 1) {
        for (plotGroup in plotGroups) {
          y = y + length(groupList[[plotGroup]]) + 1
          
          message(y)
          segments(
            x0 = 0,
            y0 = 1 - ((y - 0.5) / ncol(hic.mat.z)),
            x1 = 1,
            y1 = 1 - ((y - 0.5) / ncol(hic.mat.z)),
            lwd = 2
          )
          
          #abline(h=50.5, col="red", lwd=2)  # Draw a horizontal line
          
        }
      }
      
      
      
      
      if (middleLine) {
        segments(
          x0 = 0.5,
          y0 = 0,
          x1 = 0.5,
          y1 = 1
        )
      }
      
      
    } else{
      par(mar = c(0.5, 2, 0, 0) + 0.1)  #c(bottom, left, top, right).
      # if(subSampleGW){
      #   message('subsample to same number of rows as TetO')
      #   maxTetOrows<-max(length(TetOrow_upfw),length(TetOrow_uprv),length(TetOrow_downfw),length(TetOrow_downrv))
      #   subsetID<-floor(seq(from =1, to =length(GWrow),length.out=maxTetOrows))
      #   GWrow_subset<-GWrow[subsetID]
      # }else{
      #   GWrow_subset<-GWrow
      # }
      #
      #
      # #hicMatz_tetO<-hic.mat.z[TetOrow,c(minBin:maxBin)]
      #
      # hicMatz_tetO_upfw<-hic.mat.z[TetOrow_upfw,c(minBin:maxBin)]
      # hicMatz_tetO_uprv<-hic.mat.z[TetOrow_uprv,c(minBin:maxBin)]
      # hicMatz_tetO_dwnfw<-hic.mat.z[TetOrow_downfw,c(minBin:maxBin)]
      # hicMatz_tetO_dwnrv<-hic.mat.z[TetOrow_downrv,c(minBin:maxBin)]
      # hicMatz_GW  <-hic.mat.z[GWrow_subset,c(minBin:maxBin)]
      #
      # message('matrix1')
      #
      # hic.mat.z<-rotate(hicMatz_tetO_upfw)
      # image(hic.mat.z, col=myColor, axes=FALSE, zlim=c(0,cutoff))
      # if(middleLine){
      #   segments(x0 = 0.5, y0=0, x1 = 0.5, y1 = 1)
      # }
      # box(lwd = 2)
      #
      # message('matrix2')
      # hic.mat.z<-rotate(hicMatz_tetO_uprv)
      # image(hic.mat.z, col=myColor, axes=FALSE, zlim=c(0,cutoff))
      # if(middleLine){
      #   segments(x0 = 0.5, y0=0, x1 = 0.5, y1 = 1)
      # }
      # box(lwd = 2)
      #
      # message('matrix3')
      # hic.mat.z<-rotate(hicMatz_tetO_dwnfw)
      # image(hic.mat.z, col=myColor, axes=FALSE, zlim=c(0,cutoff))
      # if(middleLine){
      #   segments(x0 = 0.5, y0=0, x1 = 0.5, y1 = 1)
      # }
      # box(lwd = 2)
      #
      # message('matrix4')
      # hic.mat.z<-rotate(hicMatz_tetO_dwnrv)
      # image(hic.mat.z, col=myColor, axes=FALSE, zlim=c(0,cutoff))
      # if(middleLine){
      #   segments(x0 = 0.5, y0=0, x1 = 0.5, y1 = 1)
      # }
      # box(lwd = 2)
      #
      # message('matrix5-GW')
      # hic.mat.z<-rotate(hicMatz_GW)
      # image(hic.mat.z, col=myColor, axes=FALSE, zlim=c(0,cutoff))
      # if(middleLine){
      #   segments(x0 = 0.5, y0=0, x1 = 0.5, y1 = 1)
      # }
      # box(lwd = 2)
      # par(mar=c(0, 2, 0.5, 0) +0.1)  #c(bottom, left, top, right).
      #
      
    }
    
    par(mar = c(2, 2, 0, 0) + 0.1)  #c(bottom, left, top, right) The default is c(5, 4, 4, 2) + 0.1.
    #plot(NULL ,xlim=c(minBin,maxBin), ylim=c(0,0), ylab="", xlab="Distance from peak center (kb)", axes = FALSE
    plot(
      NULL ,
      xlim = c(minBin, maxBin),
      ylim = c(0, 0),
      ylab = "",
      xlab = "",
      axes = FALSE
      ,
      frame.plot = FALSE,
      cex.axis = 0.9,
      cex.lab = 1
    )
    
    axis(
      side = 1,
      at = c(minBin, (maxBin + minBin) / 2, maxBin),
      labels = c(
        paste0("-", plotZoom / 1e3, "Kb"),
        paste0("0"),
        paste0("+", plotZoom / 1e3, 'Kb')
      )
    )
    
    
    
    #colourbar
    
    
    par(mar = c(0, 2, 1, 0) + 0.1)  #c(bottom, left, top, right) The default is c(5, 4, 4, 2) + 0.1.
    # scale = (length(myColor)-1)/cutoff
    # ticks = seq(0, cutoff, len=2)
    # plot(c(0,cutoff), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
    # message(5)
    # axis(1,  round(ticks,2), las=2)
    #   for (i in 1:(length(myColor)-1)) {
    #     x = (i-1)/scale
    #     rect(x,0,x+1/scale,10, col=myColor[i], border=NA)
    #   }
    #
    
    # cutoff=0.8
    #   bgColour<-'#f0f0f0'
    #   protColour<-'#636363'
    # paletteLength <- 100 #20
    # myColor <- colorRampPalette(c(bgColour,protColour))(paletteLength)
    # myBreaks <- seq(0,cutoff,length.out=paletteLength)
    #
    #
    #
    
    colMat <-
      matrix(data = myBreaks,
             nrow = paletteLength,
             byrow = TRUE)
    if (useRaster) {
      

      image(
        colMat,
        col = myColor,
        axes = FALSE,
        zlim = c(minCutoff, cutoff),
        useRaster = TRUE
      )
    } else{
      image(
        colMat,
        col = myColor,
        axes = FALSE,
        zlim = c(minCutoff, cutoff),
        useRaster = FALSE
      )
    }
    
    axis(side = 1,
         at = c(0, 1),
         labels = c(round(minCutoff,1), round(cutoff, 1)))
    
  }
  
  
  
  if (PNG | PDF) {
    dev.off()
    
   }
}


GetTornadoCov_V4<-function(region, ChIPdf, orderCov=ChIPdf$name, orderby='FLAG_ON', zoom=5000, binWidth=50, Matrix=TRUE, makeList=FALSE, bincoverage=TRUE, prefix, bigwig='pval'){
  
  debug=FALSE
  if(debug){
    
    region=flag_opt_mrgdbyhighestPeak
    ChIPdf=test_ChIPdf
    orderCov=ChIPdf$name
    orderby='C17_FLAG_ON'
    zoom=5000
    binWidth=50
    Matrix=FALSE
    makeList=FALSE
    bincoverage=TRUE
    orderCovname=orderCov[1]
    
    
  }
  
  stopifnot(is(region, "GRanges"))
  stopifnot(is(ChIPdf, "data.frame"))
  
  infoDF<-data.frame(region=deparse(substitute(region)),
                     orderCov=paste(ChIPdf$name,collapse  = '_'),
                     orderby=ifelse(is.null(orderby),'Not ordered',orderby),
                     zoom=zoom,
                     binWidth=binWidth,
                     Matrix=Matrix,
                     makeList=makeList,
                     bincoverage=bincoverage,
                     bigwig=bigwig
  )
  
  #orderCov is now just a long string of all the ChIPdf names
  
  if(length(orderCov)>0){
    for(orderCovname in orderCov){
      if(bigwig=='pval'){
      orderChIPbw<-paste0(prefix,ChIPdf[ChIPdf$name==orderCovname,'bigwig'])
    }

    if(bigwig=='fc'){
      orderChIPbw<-paste0(prefix,ChIPdf[ChIPdf$name==orderCovname,'fc_bw'])
    }

      message('Get coverage for ',orderCovname)
      message('bw:',orderChIPbw)
      bw_score<-import(orderChIPbw, selection=BigWigSelection(region))
      gr.data.cov <- GenomicRanges::coverage(bw_score, weight="score")
      seqlevels(region) <- names(gr.data.cov)
      region<-binnedAverage(region, gr.data.cov, paste0(orderCovname,"_cov"))
    }
  }
  
  if(!is.null(orderby)){
    #Check whether coverage is already calculated for order factor otherwise get the coverage
    
    if(bigwig=='pval'){
      orderChIPbw<-paste0(prefix,ChIPdf[ChIPdf$name==orderby,'bigwig'])
    }

    if(bigwig=='fc'){
      orderChIPbw<-paste0(prefix,ChIPdf[ChIPdf$name==orderby,'fc_bw'])
    }
    

    message('Order by ',orderby)
    message('bw:',orderChIPbw)
    
    if(!paste0(orderby,"_cov") %in% colnames(mcols(region))){
      message('Get coverage for ',orderby)
      bw_score<-import(orderChIPbw, selection=BigWigSelection(region))
      gr.data.cov <- GenomicRanges::coverage(bw_score, weight="score")
      seqlevels(region) <- names(gr.data.cov)
      region<-binnedAverage(region, gr.data.cov, paste0(orderby,"_cov"))  
    }
    
    o = order(mcols(region[,paste0(orderby,'_cov')]), decreasing = TRUE)
    region<-region[o]
  }
  
  
  #Extend genomic region and tile. This will result in a GRangesList
  message('Tile regions')
  ROI<-GenomicRanges::resize(region,  width = zoom, fix = 'center')
  ROI_tiles<-GenomicRanges::tile(x=ROI, width=binWidth)
  tiles.gr <- unlist(ROI_tiles)
  
  
  if(Matrix){
    Mat_list<-list()
  }
  
  for(a in 1:nrow(ChIPdf)){

  if(bigwig=='pval'){
        message('Get coverage for ',ChIPdf$name[a],' from ', ChIPdf$bigwig[a])
        bw_score<-import(paste0(prefix,ChIPdf$bigwig[a]), selection=BigWigSelection(tiles.gr))    
    }

    if(bigwig=='fc'){
      
        message('Get coverage for ',ChIPdf$name[a],' from ', ChIPdf$fc_bw[a])
        bw_score<-import(paste0(prefix,ChIPdf$fc_bw[a]), selection=BigWigSelection(tiles.gr))    

    }



    
    if(bincoverage){
      message('Calculating averge bin coverage')
      gr.data.cov <- GenomicRanges::coverage(bw_score, weight="score")
      seqlevels(tiles.gr) <- names(gr.data.cov)
      tiles.gr<-binnedAverage(tiles.gr, gr.data.cov, paste0(ChIPdf$name[a],"_cov"))
    }else{
      message('Calculating average bin score')
      hits <- findOverlaps(tiles.gr, bw_score)
      agg <- aggregate(bw_score, hits, score=mean(score))
      mcols(tiles.gr)[paste0(ChIPdf$name[a],"_cov")] <- agg$score
    }
    if(Matrix){
      Mat_list[[ChIPdf$prot[a]]][[ChIPdf$condition[a]]]<-matrix(data = mcols(tiles.gr)[,paste0(ChIPdf$name[a],'_cov')], nrow=length(region), byrow = TRUE)
    }
    
  }
  
  
  if(Matrix){
    return(list(peaks=region, Mat=Mat_list, ChIPdf=ChIPdf, infoDF=infoDF))
  }else{
    
    if(makeList){
      tiles.gr<-relist(tiles.gr, ROI_tiles)
      return(list(peaks=region, tiles=tiles.gr, ChIPdf=ChIPdf, infoDF=infoDF))
    }else{
      return(list(peaks=region, tiles=tiles.gr, ChIPdf=ChIPdf,infoDF=infoDF))
    }
    
  }
  
  
  
  
}





TACL_overlay4CfromLocalNormalized <- function(normalized_4C, exp1, exp2, plotZoom, yMax = 2500, name1 = NULL, name2 = NULL, vLine = NULL) {
  gr <- subsetByOverlaps(normalized_4C, plotZoom)


  DF <- data.frame(pos = gr$pos, V4C1 = mcols(gr)[[paste0("norm4C_", exp1)]], V4C2 = mcols(gr)[[paste0("norm4C_", exp2)]])
  DF$V4Cmin <- apply(DF[, 2:3], 1, min)
  DF$V4Cmax <- apply(DF[, 2:3], 1, max)
  DF$colors <- ifelse(DF$V4C1 > DF$V4C2, "forestgreen", "orange")

  # Plotting
  plot(
    x = DF$pos, y = DF$V4Cmax, type = "h", col = DF$colors, frame.plot = FALSE,
    axes = FALSE,
    xlim = c(start(plotZoom), end(plotZoom)), ylim = c(0, yMax), ann = FALSE
  )
  points(x = DF$pos, y = DF$V4Cmin, type = "h", col = "lightgray", ylim = c(0, yMax))


    
  if (!is.null(vLine)) {
      abline(v = vLine, col = "red", lty = 2)
  }

  axis(2, seq(0, yMax, yMax), las = 2)
  mtext("4C", side = 2, line = 3.5, at = yMax / 2, adj = 1, cex = 0.7, las = 1)

  if (is.null(name1)) {
    name1 <- exp1
  }
  if (is.null(name2)) {
    name2 <- exp2
  }

  legend("topright",
    legend = c(name1, name2),
    fill = c("forestgreen", "orange"),
    bty = "n"
  )
}


plotCTCFmotif <- function(peaks, zoom, y1 = -5, y2 = -10, strandCol = "FIMO", Mb = TRUE, scaleTriangle = 2.5e-08) {
  # Check if the peaks object is of class GRanges
  if (!inherits(peaks, "GRanges")) {
    stop("ERROR: peaks file is not a GRanges object")
  }

  # Check if the specified strand column exists in the metadata columns of peaks
  if (!strandCol %in% colnames(mcols(peaks))) {
    stop(paste("ERROR: peaks file does not contain column", strandCol))
  }

  x.wid <- scaleTriangle * width(zoom)

  # x.wid <- 2.5e-16 * width(zoom)

  if (Mb) {
    Xdiv <- 1e6
  } else {
    Xdiv <- 1
  }

  MOI <- subsetByOverlaps(peaks, zoom)
  MOI_withMotif <- MOI[as.vector(!is.na(mcols(MOI)[paste0(strandCol)]))]


  # plot positive strands
  positive_strand_peaks <- MOI_withMotif[mcols(MOI_withMotif)[[strandCol]] == "+"]

  if (length(positive_strand_peaks) > 0) {
    triangle <- lapply(start(positive_strand_peaks) / Xdiv, function(x) {
      list(
        xx = c(x, x, x + x.wid),
        yy = c(y1, y2, (y1 + y2) / 2)
      )
    })
    lapply(triangle, function(x) polygon(x$xx, x$yy, col = "red", border = "red"))
  }

  # plot negative strands
  negative_strand_peaks <- MOI_withMotif[mcols(MOI_withMotif)[[strandCol]] == "-"]
  if (length(negative_strand_peaks) > 0) {
    triangle <- lapply(start(negative_strand_peaks) / Xdiv, function(x) {
      list(
        xx = c(x, x, x - x.wid),
        yy = c(y1, y2, (y1 + y2) / 2)
      )
    })
    lapply(triangle, function(x) polygon(x$xx, x$yy, col = "blue", border = "blue"))
  }

  # plot unorientated motifs as grey rectangles
  both_strand_peaks <- MOI_withMotif[mcols(MOI_withMotif)[[strandCol]] == "both"]
  if (length(both_strand_peaks) > 0) {
    rect(start(both_strand_peaks) / Xdiv, y1, end(both_strand_peaks) / Xdiv, y2, col = "lightgrey", border = "lightgrey")
  }

  # plot strands without a motif as black rectangles
  noMotif_peaks <- MOI[is.na(mcols(MOI)[[strandCol]])]
  if (length(noMotif_peaks) > 0) {
    rect(start(noMotif_peaks) / Xdiv, y1, end(noMotif_peaks) / Xdiv, y2, col = "black", border = "black")
  }
}


plot_bw <- function(bw, zoom, name = "", expName = "", yMax = 100, colour = "grey", peaks = NULL, scaleFactor = 1, BWlwd=1.5) {
  bw <- rtracklayer::import(bw, which = zoom)
  bw$col <- colour

  if (!is.null(peaks)) {
    IDX <- findOverlaps(bw, peaks)@from
    if (length(IDX) > 0) {
      bw[IDX]$col <- "red"
    }
  }

  plot(start(bw) / 1e6, bw$score/scaleFactor,
    ylim = c(0, yMax),
    xlim = c(start(zoom) / 1e6, end(zoom) / 1e6),
    xlab = "", ylab = "", type = "h",
    frame.plot = FALSE, cex.axis = 0.9, cex.lab = 1,
    cex.main = 0.9,
    col = bw$col,
    axes = FALSE,
    ann = FALSE, 
    lwd=BWlwd
  )

  axis(2, c(0, yMax), las = 2)
  mtext(name, side = 2, line = 3.5, at = yMax / 2, adj = 1, cex = 0.7, las = 1)

  legend("topright", legend = expName, fill = colour, bty = "n")
}

plot_bw_TACL_expname <- function(expname, zoom, ChIPinfo, name = NULL, yMax = 100, norm='none', colour = NULL, BWlwd=2) {
  if (is.null(name)) {
    name <- ChIPinfo[ChIPinfo$name == expname, "prot"]
  }

  if (is.null(colour)) {
    colour <- protCols[protCols$protein == name, "color"]
  }

  if(norm=='none'){
    scaleFactor=1
  }

  if(norm=='meanBG'){
    scaleFactor=ChIPinfo[ChIPinfo$name == expname, "meanBG"]
  }

  if(norm=='outerBG'){
    scaleFactor=ChIPinfo[ChIPinfo$name == expname, "outerBG"]
  }

  plot_bw(file.path("/storage/shared/", ChIPinfo[ChIPinfo$name == expname, "bigwig"]), zoom, name = name, expName = expname, yMax = yMax, colour = colour, scaleFactor = scaleFactor, BWlwd=BWlwd)
}

