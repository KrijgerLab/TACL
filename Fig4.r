
#fIGURE 4
#4C/D in degron backgrounds ------------------------------------------------



# Load packages -----------------------------------------------------------
library("GenomicRanges")
require("caTools")

library(foreach)
library(doParallel)


# functions ----------------------------------------------------------------

prefix <- ""
source(file.path(prefix, "./TACL_functions_NG2025.r"))


# features ----------------------------------------------------------------

folder_4C <- "./normalized_4c"   #Github page
folder_peaks <- "./peaks"  #Github page

outF <- "./out/Fig4"
if (!dir.exists(outF)) {
  dir.create(outF, recursive = TRUE)
}

# functions ---------------------------------------------------------------

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


TACL_overlay4CfromLocalNormalized<-function(normalized_4C, exp1, exp2, plotZoom, yMax=2500, name1=NULL, name2=NULL){


  gr<-subsetByOverlaps (normalized_4C, plotZoom)
 
  
  DF <- data.frame( pos = gr$pos, V4C1 = mcols(gr)[[paste0("norm4C_", exp1)]], V4C2 = mcols(gr)[[paste0("norm4C_", exp2)]] )
  DF$V4Cmin <- apply( DF[,2:3], 1, min )
  DF$V4Cmax <- apply( DF[,2:3], 1, max )
  DF$colors <- ifelse( DF$V4C1 > DF$V4C2, "forestgreen", "orange" )
  
  #Plotting
  plot( x=DF$pos, y=DF$V4Cmax, type='h', col=DF$colors, frame.plot=FALSE
        ,axes = FALSE
        ,xlim=c(start(plotZoom),end(plotZoom)), ylim=c(0,yMax), ann=FALSE)
  points( x=DF$pos, y=DF$V4Cmin, type='h', col='lightgray', ylim=c(0,yMax) )
  axis(2, seq(0,yMax,yMax),las=2)


  if(is.null(name1)){
    name1<-exp1
  }
  if(is.null(name2)){
    name2<-exp2
  }

  legend("topright"
         ,legend=c(name1,name2)
         ,fill=c("forestgreen", "orange")
         ,bty ="n")



}


# #First extract local captures

# list4C <- list()

# list4C$WT <- local4C(readRDS("/storage/shared/TACL/PK/4C/VER10253/VER10253_TACL_MR_4C_eHAP1_C17_WT_TetO_rep1.rds"), ROI = C17TetO, normZoom = 2e7)
# list4C$on <- local4C(readRDS("/storage/shared/TACL/GEO/rds/4Cseq_TMAU2.rds"), ROI = C17TetO, normZoom = 2e7)
# list4C$off <- local4C(readRDS("/storage/shared/TACL/GEO/rds/4Cseq_TMAU2_Dox1h.rds"), ROI = C17TetO, normZoom = 2e7)

# list4C$WAPL_on_IAA <- local4C(readRDS("/storage/shared/TACL/data/4C/VER9308//VER9308_TACL_MR_4C_eHAP1_C17_WAPLAID_C8_TMAU2TIR1_IAA3h_4C_rep1.rds"), ROI = C17TetO, normZoom = 2e7)
# list4C$WAPL_off_IAA <- local4C(readRDS("/storage/shared/TACL/data/4C/VER9308//VER9308_TACL_MR_4C_eHAP1_C17_WAPLAID_C8_TMAU2TIR1_IAA3h_Dox3h_4C_rep1.rds"), ROI = C17TetO, normZoom = 2e7)
# list4C$WAPL_IAA <- local4C(readRDS("/storage/shared/TACL/PK/4C/VER10253/VER10253_TACL_MR_4C_eHAP1_C17_WAPLAID_C11_IAA3h_TetO_rep1.rds"), ROI = C17TetO, normZoom = 2e7)

# list4C$CTCF_on_IAA <- local4C(readRDS("/storage/shared/TACL/data/4C/VER9308//VER9308_TACL_MR_4C_eHAP1_C17_CTCFAID_C22_TMAU2TIR1_IAA3h_4C_rep1.rds"), ROI = C17TetO, normZoom = 2e7)
# list4C$CTCF_off_IAA <- local4C(readRDS("/storage/shared/TACL/data/4C/VER9308//VER9308_TACL_MR_4C_eHAP1_C17_CTCFAID_C22_TMAU2TIR1_IAA3h_Dox3h_4C_rep1.rds"), ROI = C17TetO, normZoom = 2e7)
# list4C$CTCF_IAA <- local4C(readRDS("/storage/shared/TACL/PK/4C/VER10253/VER10253_TACL_MR_4C_eHAP1_C17_CTCFAID_C22_IAA3h_TetO_rep1.rds"), ROI = C17TetO, normZoom = 2e7)

# list4C$STAG2_on_IAA <- local4C(readRDS("/storage/shared/TACL/data/4C/VER9308//VER9308_TACL_MR_4C_eHAP1_C17_STAG2AID_C2_TMAU2TIR1_IAA1h_4C_rep1.rds"), ROI = C17TetO, normZoom = 2e7)
# list4C$STAG2_IAA <- local4C(readRDS("/storage/shared/TACL/PK/4C/VER10253/VER10253_TACL_MR_4C_eHAP1_C17_STAG2AID_C2_IAA1h_TetO_rep1.rds"), ROI = C17TetO, normZoom = 2e7)

# list4C$PDS5A_on_IAA <- local4C(readRDS("/storage/shared/TACL/data/4C/VER9308//VER9308_TACL_MR_4C_eHAP1_C17_PDS5AAID_B10_TMAU2TIR1_IAA3h_4C_rep1.rds"), ROI = C17TetO, normZoom = 2e7)
# list4C$PDS5A_IAA <- local4C(readRDS("/storage/shared/TACL/PK/4C/VER10253/VER10253_TACL_MR_4C_eHAP1_C17_PDS5AAID_B10_IAA3h_TetO_rep1.rds"), ROI = C17TetO, normZoom = 2e7)

# # Normalize 4C data -------------------------------------------------------
# normalized_4C <- normalize4C(list4C, listNames = c("WT", "on", "off", "WAPL_on_IAA", "WAPL_off_IAA", "WAPL_IAA", "CTCF_on_IAA", "CTCF_off_IAA", "CTCF_IAA", "STAG2_on_IAA", "STAG2_IAA", "PDS5A_on_IAA", "PDS5A_IAA"), ROI = C17TetO, normZoom = 2e7, wSize = 21)
   

# saveRDS(normalized_4C, "/storage/shared/TACL/PK/4C/TACL_C17_degrons_including_VER10253_normalized_4C.rds")




normalized_4C_VER10253<- readRDS("./normalized_4C/TACL_C17_degrons_including_VER10253_normalized_4C.rds")


# listNames = c("WT", "on", "off", "WAPL_on_IAA", "WAPL_off_IAA", "WAPL_IAA", "CTCF_on_IAA", "CTCF_off_IAA", "CTCF_IAA", "STAG2_on_IAA", "STAG2_IAA", "PDS5A_on_IAA", "PDS5A_IAA"), ROI = C17TetO, normZoom = 2e7, wSize = 21)




plot4Cs_degron <- function(a, C17TetO, ChIPinfo, domain_HMM_GR, outF, protCols, normalized_4C, yMax4C=750, zoomWidth=8e6, PNG = FALSE, PDF = FALSE) {
    message(a)

    if (!dir.exists(outF)) {
        dir.create(outF, recursive = TRUE)
    }

    ROI <- C17TetO[a]
    smallDomain <- subsetByOverlaps(domain_HMM_GR, ROI)

    zoom <- resize(smallDomain, width = zoomWidth, fix = "center")

        #zoom <- resize(smallDomain, width = width(smallDomain) + flankingWidth, fix = "center")


    height_vector <- c(rep(lcm(3), 9), lcm(2), rep(lcm(0.5), 4))
    height_layout <- sum(as.numeric(gsub("cm", "", height_vector)))
    if (PNG) {
        png(filename = paste0(outF, "/Degron4C_TetO_", a, "_", seqnames(ROI), "_HMMdomain_localNormalized.png"), width = 25, height = height_layout + 3, units = "cm", res = 100)
    }
    if (PDF) {
        pdf(file = paste0(outF, "/Degron4C_TetO_", a, "_", seqnames(ROI), "_HMMdomain_localNormalized.pdf"), width = 25 / 2.54, height = (height_layout + 3) / 2.54)

    }

    mat.row <- length(height_vector)
    layout_matrix <- matrix(1:mat.row, nrow = mat.row, ncol = 1, byrow = TRUE)
    layout(
        mat = layout_matrix,
        height = height_vector,
        width = c(rep(lcm(20), mat.row))
    )

    par(oma = c(0, 0, 4, 0))
    par(mar = c(0, 6, 0, 2) + 0.1)


    # 4C plots (9x)


    #TACL_overlay4CfromLocalNormalized(normalized_4C[[a]], exp1 = "off", exp2 = "WT", plotZoom = zoom, yMax = yMax4C, name1 = "TACL OFF", name2 = "WT")
    TACL_overlay4CfromLocalNormalized(normalized_4C[[a]], exp1 = "off", exp2 = "on", plotZoom = zoom, yMax = yMax4C, name1 = "TACL OFF", name2 = "TACL ON")
    
    #WAPL
    TACL_overlay4CfromLocalNormalized(normalized_4C[[a]], exp1 = "off", exp2 = "WAPL_IAA", plotZoom = zoom, yMax = yMax4C, name1 = "off", name2 = "WAPL depleted")
    TACL_overlay4CfromLocalNormalized(normalized_4C[[a]], exp1 = "WAPL_IAA", exp2 = "WAPL_on_IAA", plotZoom = zoom, yMax = yMax4C, name1 = "WAPL depleted", name2 = "TACL-ON; WAPL depleted")
    
    #CTCF
    TACL_overlay4CfromLocalNormalized(normalized_4C[[a]],
        exp1 = "off",
        exp2 = "CTCF_IAA",
        plotZoom = zoom,
        yMax = yMax4C,
        name1 = "off",
        name2 = "CTCF depleted")

    TACL_overlay4CfromLocalNormalized(normalized_4C[[a]],
        exp1 = "CTCF_IAA",
        exp2 = "CTCF_on_IAA",
        plotZoom = zoom,
        yMax = yMax4C,
        name1 = "CTCF depleted",
        name2 = "TACL-ON; CTCF depleted")

    #PDS5A
    TACL_overlay4CfromLocalNormalized(normalized_4C[[a]],
        exp1 = "off",
        exp2 = "PDS5A_IAA",
        plotZoom = zoom,
        yMax = yMax4C,
        name1 = "off",
        name2 = "PDS5A depleted")

    TACL_overlay4CfromLocalNormalized(normalized_4C[[a]],
        exp1 = "PDS5A_IAA",
        exp2 = "PDS5A_on_IAA",
        plotZoom = zoom,
        yMax = yMax4C,
        name1 = "PDS5A depleted",
        name2 = "TACL-ON; PDS5A depleted")

    
    #STAG2
    TACL_overlay4CfromLocalNormalized(normalized_4C[[a]],
        exp1 = "off",
        exp2 = "STAG2_IAA",
        plotZoom = zoom,
        yMax = yMax4C,
        name1 = "off",
        name2 = "STAG2 depleted")

    TACL_overlay4CfromLocalNormalized(normalized_4C[[a]],
        exp1 = "STAG2_IAA",
        exp2 = "STAG2_on_IAA",
        plotZoom = zoom,
        yMax = yMax4C,
        name1 = "STAG2 depleted",
        name2 = "TACL-ON; STAG2 depleted")


   
    # CTCF
    plot_bw_TACL_expname(expname = "C17_TACL-ON_CTCF_rep1", zoom, ChIPinfo, yMax = 3, norm = "meanBG", BWlwd = 1)
 

    # CTCF motif
    CTCF_zoom <- peak_CTCF_C17_TMAU2[findOverlaps(peak_CTCF_C17_TMAU2, zoom)@from]
    plot(NA, type = "n", xlim = c(start(zoom) / 1e6, end(zoom) / 1e6), ylim = c(0, 1), axes = FALSE, ann = FALSE)
    if (length(CTCF_zoom) > 0) {
        plotCTCFmotif(peaks = CTCF_zoom, zoom, y1 = 0, y2 = 1, strandCol = "FIMO", Mb = TRUE, scaleTriangle = 1e-08)
    }
    mtext("CTCF", side = 2, line = 3.5, adj = 1, cex = 0.7, las = 1)


 
    # TetO location
    TetO_zoom <- C17TetO[findOverlaps(C17TetO, zoom)@from]
    plot(NA, type = "n", xlim = c(start(zoom), end(zoom)), ylim = c(0, 1), axes = FALSE, ann = FALSE)
    if (length(TetO_zoom) > 0) {
        rect(xleft = start(TetO_zoom), xright = end(TetO_zoom), ybottom = 0, ytop = 1, col = "darkgreen", border = "darkgreen")
    }
    mtext("TetO", side = 2, line = 3.5, adj = 1, cex = 0.7, las = 1)

    # HMM domain
    plot(NA, type = "n", xlim = c(start(zoom), end(zoom)), ylim = c(0, 1), axes = FALSE, ann = FALSE)
    if (length(smallDomain) > 0) {
        rect(xleft = start(smallDomain), xright = end(smallDomain), ybottom = 0, ytop = 1, col = "darkred", border = "darkred")
    }
    mtext("HMM domain", side = 2, line = 3.5, adj = 1, cex = 0.7, las = 1)


     # X-axis
    par(mar = c(0.5, 6, 0, 2) + 0.1)
    plot(NA, type = "n", xlim = c(start(zoom), end(zoom)), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)
    axis(1, at = seq(from = start(zoom), to = end(zoom), by = (end(zoom) - start(zoom)) / 5), las = 1)


    # Title
    mtext(ROI, outer = TRUE, cex = 1, line = -0.5)
    if (PNG | PDF) {
        dev.off()
    }
}


# #test
# plot4Cs_degron(a = 21, C17TetO, ChIPinfo = ChIPdf, domain_HMM_GR, outF = "/storage/shared/TACL/PK/manuscript/revision/fig/STAG/4C/STAG2_peaks_domain_V4/", protCols = protCols, normalized_4C=normalized_4C_VER10253, yMax4C=750, zoomWidth=12e6, PNG = FALSE, PDF = FALSE)



#zoom
library(foreach)
library(doParallel)

numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)


foreach(a = seq_along(C17TetO), .packages = c("GenomicRanges", "graphics", "grDevices", "grid")) %dopar% {
    plot4Cs_degron(a , C17TetO, ChIPinfo = ChIPdf, domain_HMM_GR, outF = outF, protCols = protCols, normalized_4C=normalized_4C_VER10253, yMax4C=750, zoomWidth=8e6, PNG = FALSE, PDF = TRUE)
}

# Stop the cluster
stopCluster(cl)
