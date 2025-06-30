# R functions to normalize and plot 4C data around TetO sites in the C17 cell line.

# Load packages -----------------------------------------------------------
library("GenomicRanges")
require("caTools")

library(foreach)
library(doParallel)


# functions ----------------------------------------------------------------

# Add all functions to source file

prefix <- ""
source(file.path(prefix, "./TACL_functions_NG2025.r"))


# features ----------------------------------------------------------------

folder_4C <- "./normalized_4c"   #Github page
folder_peaks <- "./peaks"  #Github page

outF <- "./out/Fig1"
if (!dir.exists(outF)) {
  dir.create(outF, recursive = TRUE)
}

# Fig 1D: 4C example plot with CTCF and TetO sites ----------------------------------------------------------------


# Make a reduced 4C coverage file with only coverage within 10MB of the TetO site



#output file can be downloaded from the github page normalized_4C folder
  rds_mch_2 <- readRDS(file.path(folder_4C, "Cherry.rds"))
  rds_TACL_2 <- readRDS(file.path(folder_4C, "TACL-ON.rds"))
  rds_dox_2 <- readRDS(file.path(folder_4C, "TACL-OFF.rds"))

  list4C <- list()
  list4C$mch_2 <- local4C(rds_mch_2, ROI = C17TetO, normZoom = 2e7)
  list4C$on_2 <- local4C(rds_TACL_2, ROI = C17TetO, normZoom = 2e7)
  list4C$off_2 <- local4C(rds_dox_2, ROI = C17TetO, normZoom = 2e7)

  # Normalize the 4C data
  normalized_4C <- normalize4C(list4C, listNames = c("on_2", "off_2", "mch_2"), ROI = C17TetO, normZoom = 2e7, wSize = 21)
  saveRDS(normalized_4C, file.path(folder_4C, "TACL_rep2_normalized_4C.rds"))


# Plot delta 4C plots, with CTCF ChIPseq track

normalized_4C <- readRDS(file.path(folder_4C, "TACL_rep2_normalized_4C.rds"))
peak_CTCF_C17_TMAU2 <- readRDS(file.path(folder_4C, "C17_TACL_CTCF_hg38_peaks_withMotif_TetO.rds"))
peak_CTCF_C17_Tmcherry <- readRDS(file.path(folder_4C, "C17_mCherry_CTCF_hg38_peaks_withMotif_TetO.rds"))




plot4Cs_fig1D <- function(a, C17TetO, ChIPinfo, domain_HMM_GR, zoomWidth, scaleTriangle, outF, pdf = FALSE, png = FALSE) {
  message(a)

  if (!dir.exists(outF)) {
    dir.create(outF, recursive = TRUE)
  }

  ROI <- C17TetO[a]
  zoom <- resize(ROI, width = zoomWidth, fix = "center")

  height_vector <- c(rep(lcm(4), 2), rep(lcm(2), 1), lcm(0.5), lcm(2), rep(lcm(0.5), 3))
  height_layout <- sum(as.numeric(gsub("cm", "", height_vector)))

  if (png) {
    png(filename = paste0(outF, "/off_vs_degrons_TetO_", a, "_", seqnames(ROI), "_", zoomWidth / 1e6, "MB_localNormalized.png"), width = 25, height = height_layout + 3, units = "cm", res = 100)
  }

  if (pdf) {
    pdf(file = paste0(outF, "/off_vs_degrons_TetO_", a, "_", seqnames(ROI), "_", zoomWidth / 1e6, "MB_localNormalized.pdf"), width = 25 / 2.54, height = (height_layout + 3) / 2.54)
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


  # layout.show(mat.row)

  # 4C plots (2x)
  yMax4C <- 750
  TACL_overlay4CfromLocalNormalized(normalized_4C[[a]], exp1 = "mch_2", exp2 = "on_2", plotZoom = zoom, yMax = yMax4C, name1 = "TACL-Cherry", name2 = "TACL-ON")
  TACL_overlay4CfromLocalNormalized(normalized_4C[[a]], exp1 = "off_2", exp2 = "on_2", plotZoom = zoom, yMax = yMax4C, name1 = "TACL-OFF", name2 = "TACL-ON")


  # CTCF
  plot_bw(file.path("/storage/shared/", ChIPinfo[ChIPinfo$name == "C17_mCherry_CTCF_rep1", "bigwig"]), zoom, name = "CTCF", expName = "C17_mCherry_CTCF_rep1", yMax = 300, colour = protCols[protCols$protein == "CTCF", "color"])

  # CTCF motif

  CTCF_zoom <- peak_CTCF_C17_Tmcherry[findOverlaps(peak_CTCF_C17_Tmcherry, zoom)@from]
  plot(NA, type = "n", xlim = c(start(zoom) / 1e6, end(zoom) / 1e6), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  if (length(CTCF_zoom) > 0) {
    plotCTCFmotif(peaks = CTCF_zoom, zoom, y1 = 0, y2 = 1, strandCol = "FIMO", Mb = TRUE, scaleTriangle = scaleTriangle)
  }
  mtext("CTCF", side = 2, line = 3.5, adj = 1, cex = 0.7, las = 1)


  plot_bw(file.path("/storage/shared/", ChIPinfo[ChIPinfo$name == "C17_TACL-ON_CTCF_rep1", "bigwig"]), zoom, name = "CTCF", expName = "C17_TACL-ON_CTCF_rep1", yMax = 300, colour = protCols[protCols$protein == "CTCF", "color"])

  # CTCF motif
  CTCF_zoom <- peak_CTCF_C17_TMAU2[findOverlaps(peak_CTCF_C17_TMAU2, zoom)@from]
  plot(NA, type = "n", xlim = c(start(zoom) / 1e6, end(zoom) / 1e6), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  if (length(CTCF_zoom) > 0) {
    plotCTCFmotif(peaks = CTCF_zoom, zoom, y1 = 0, y2 = 1, strandCol = "FIMO", Mb = TRUE, scaleTriangle = scaleTriangle)
  }
  mtext("CTCF", side = 2, line = 3.5, adj = 1, cex = 0.7, las = 1)



  # TetO location
  TetO_zoom <- C17TetO[findOverlaps(C17TetO, zoom)@from]
  plot(NA, type = "n", xlim = c(start(zoom), end(zoom)), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  if (length(TetO_zoom) > 0) {
    rect(xleft = start(TetO_zoom), xright = end(TetO_zoom), ybottom = 0, ytop = 1, col = "darkgreen", border = "darkgreen")
  }
  mtext("TetO", side = 2, line = 3.5, adj = 1, cex = 0.7, las = 1)

  # X-axis
  par(mar = c(0.5, 6, 0, 2) + 0.1)
  plot(NA, type = "n", xlim = c(start(zoom), end(zoom)), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)
  axis(1, at = seq(from = start(zoom), to = end(zoom), by = (end(zoom) - start(zoom)) / 5), las = 1)

  # Title
  mtext(ROI, outer = TRUE, cex = 1, line = -0.5)

  if (pdf | png) {
    dev.off()
  }
}

#test plot
# plot4Cs_fig1D(a = 3, C17TetO, ChIPinfo = ChIPdf, domain_HMM_GR, zoomWidth = 6e6, outF = "/storage/shared/TACL/PK/manuscript/revision/fig/1D/", scaleTriangle = 1e-08)


numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)


foreach(a = seq_along(C17TetO), .packages = c("GenomicRanges", "graphics", "grDevices", "grid")) %dopar% {
  plot4Cs_fig1D(a, C17TetO, ChIPinfo = ChIPdf, domain_HMM_GR, zoomWidth = 4e6, outF = paste0(outF, "/1D/"), scaleTriangle = 1e-08, png = TRUE)
}

# Stop the cluster
stopCluster(cl)



























