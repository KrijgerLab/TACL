# Figure 2

# Load packages -----------------------------------------------------------
library("GenomicRanges")
require("caTools")

library(foreach)
library(doParallel)



# functions ----------------------------------------------------------------

# Add all functions to source file


prefix <- ""
source(file.path(prefix, "./TACL_functions_NG2025.r"))
source(file.path(prefix, "./TACL_features_NG2025.r"))


outF <- "./out/Fig2"
if (!dir.exists(outF)) {
  dir.create(outF, recursive = TRUE)
}

covF <- "./cov"
if (!dir.exists(covF)) {
  dir.create(covF, recursive = TRUE)
}

# 2A FLAG coverage ON vs OFF --------------------------------
library(ggExtra)


#normalized flag counts at flag peaks
normcounts <- read.table("./input/flag_SE_normcounts.txt", header = TRUE, sep = "\t")  #github page input folder

#flag peaks
flag_opt_df <- as.data.frame(flag_opt_mrgdbyhighestPeak)[, c(1:3, 8, 12)] #seqnames   start     end signalValue ID
merged_df <- merge(normcounts[, -c(5, 6)], flag_opt_df, by = c("seqnames", "start", "end")) #merge by coordinates to be sure.


# Convert merged data frame back to GRanges
merged_GR <- makeGRangesFromDataFrame(merged_df,
    keep.extra.columns = TRUE,
    seqnames.field = c("seqnames"),
    start.field = c("start"),
    end.field = c("end")
)

merged_GR<-sort(merged_GR)


  #dFLAG <- merged_GR[merged_GR$M > 1 & (merged_GR$on_avg + 1) > 2^4.5, ]

  #dFLAG peaks within domains. 421 have a signal value > 35
  # table(dFLAG$signalValue > 35, dFLAG$domain>0)

  #         FALSE TRUE
  #   FALSE    11   33
  #   TRUE     39  421





# Define a pseudocount to avoid log(0) issues
pseudocount <- 1

# Create a data frame from the GRanges object for easier plotting with ggplot2
df <- as.data.frame(merged_GR)
df$domain_status <- ifelse(df$domain > 0, "Domain", "Non-domain")

# Define color palette
colors <- c("Non-domain" = "#0072B2", "Domain" = "#D55E00")

# Create a ggplot
p <- ggplot(df, aes(x = log2(on_avg + pseudocount), y = log2(off_rep2 + pseudocount))) +
  geom_point(aes(color = domain_status, shape = domain_status), size = 3, alpha = 0.7) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = c("Non-domain" = 16, "Domain" = 17)) +
  labs(
    title = "Log2 On vs Off Rep2 Signal",
    x = "log2(On Average)",
    y = "log2(Off Rep2)",
    color = "Domain Status",
    shape = "Domain Status"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),
    legend.position = "top",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) +
  xlim(0, 10) +
  ylim(0, 10)

# Add marginal density plots
p <- ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)

ggsave(file.path(outF, "fig2a.pdf"), p, width = 10, height = 10)




# 2B - Example coverage plots --------------------------------
normalized_4C_rep2 <- readRDS(file.path(folder_4C, "TACL_rep2_normalized_4C.rds"))


plot4Cs_fig2B_reps <- function(a, C17TetO, ChIPinfo, domain_HMM_GR, zoomWidth, scaleTriangle,  yMax4C=1000, BWnorm='meanBG', BWlwd =2, outF, pdf = FALSE, png = FALSE, normalized_4C_rep2=normalized_4C_rep2) {
  message(a)

  if (png | pdf) {
    if (!dir.exists(outF)) {
      dir.create(outF, recursive = TRUE)
    }
  }

  ROI <- C17TetO[a]
  zoom <- resize(ROI, width = zoomWidth, fix = "center")

  height_vector <- c(rep(lcm(3), 1), rep(lcm(2), 17), rep(lcm(0.5), 4))
  height_layout <- sum(as.numeric(gsub("cm", "", height_vector)))

  if (png) {
    png(filename = paste0(outF, "/2B_off_on_mch_TetO_", a, "_", seqnames(ROI), "_", zoomWidth / 1e6, "MB_localNormalized.png"), width = 25, height = height_layout + 3, units = "cm", res = 100)
  }

  if (pdf) {
    pdf(file = paste0(outF, "/2B_off_on_mch_TetO_", a, "_", seqnames(ROI), "_", zoomWidth / 1e6, "MB_localNormalized.pdf"), width = 25 / 2.54, height = (height_layout + 3) / 2.54)
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






  message("4C plots")
  TACL_overlay4CfromLocalNormalized(normalized_4C_rep2[[a]], exp1 = "off_2", exp2 = "on_2", plotZoom = zoom, yMax = yMax4C, name1 = "TACL-OFF", name2 = "TACL-ON")

  yMax <- 125
  if (BWnorm == 'meanBG') {
    yMax <- 5
  }
  plot_bw_TACL_expname(expname = "C17_TACL-ON_FLAG_rep1_rep2", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_TACL-OFF_FLAG_rep2", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  
  yMax <- 125
  if (BWnorm == 'meanBG') {
    yMax <- 5
  }
  plot_bw_TACL_expname(expname = "C17_TACL-ON_MAU2_rep1_rep2", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_TACL-OFF_MAU2_rep2", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_mCherry_MAU2_rep1_rep2", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)


  yMax <- 125
  if (BWnorm == 'meanBG') {
    yMax <- 20
  }
  plot_bw_TACL_expname(expname = "C17_TACL-ON_NIPBL_rep2", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_TACL-OFF_NIPBL_rep2", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_mCherry_NIPBL_rep1_rep2", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  
  yMax <- 250
  if (BWnorm == 'meanBG') {
    yMax <- 5
  }
  
  plot_bw_TACL_expname(expname = "C17_TACL-ON_RAD21_rep1_rep3", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_TACL-OFF_RAD21_rep3", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_mCherry_RAD21_rep1_rep3", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)

  yMax <- 250
  if (BWnorm == 'meanBG') {
    yMax <- 5
  }
  plot_bw_TACL_expname(expname = "C17_TACL-ON_SMC1_rep2_rep3", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_TACL-OFF_SMC1_rep2_rep3", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
    plot_bw_TACL_expname(expname = "C17_mCherry_SMC1_rep2_rep3", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  
  yMax <- 250
  if (BWnorm == 'meanBG') {
    yMax <- 3
  }
  plot_bw_TACL_expname(expname = "C17_TACL-ON_CTCF_rep1", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_TACL-OFF_CTCF_rep1", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)
  plot_bw_TACL_expname(expname = "C17_mCherry_CTCF_rep1", zoom, ChIPinfo, yMax = yMax, norm = BWnorm, BWlwd = BWlwd)

  CTCF_zoom <- peak_CTCF_C17_TMAU2[findOverlaps(peak_CTCF_C17_TMAU2, zoom)@from]
  plot(NA, type = "n", xlim = c(start(zoom) / 1e6, end(zoom) / 1e6), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  if (length(CTCF_zoom) > 0) {
    plotCTCFmotif(peaks = CTCF_zoom, zoom, y1 = 0, y2 = 1, strandCol = "FIMO", Mb = TRUE, scaleTriangle = scaleTriangle)
  }
  mtext("CTCF", side = 2, line = 3.5, adj = 1, cex = 0.7, las = 1)

  TetO_zoom <- C17TetO[findOverlaps(C17TetO, zoom)@from]
  plot(NA, type = "n", xlim = c(start(zoom), end(zoom)), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  if (length(TetO_zoom) > 0) {
    rect(xleft = start(TetO_zoom), xright = end(TetO_zoom), ybottom = 0, ytop = 1, col = "darkgreen", border = "darkgreen") 
  }
  mtext("TetO", side = 2, line = 3.5, adj = 1, cex = 0.7, las = 1)

  domain_HMM_zoom <- domain_HMM_GR[findOverlaps(domain_HMM_GR, zoom)@from]
  plot(NA, type = "n", xlim = c(start(zoom), end(zoom)), ylim = c(0, 1), axes = FALSE, ann = FALSE)
  if (length(domain_HMM_zoom) > 0) {
    rect(xleft = start(domain_HMM_zoom), xright = end(domain_HMM_zoom), ybottom = 0, ytop = 1, col = "darkblue", border = "darkblue")
  }
  mtext("TACL domain", side = 2, line = 3.5, adj = 1, cex = 0.7, las = 1)

  par(mar = c(0.5, 6, 0, 2) + 0.1)
  plot(NA, type = "n", xlim = c(start(zoom), end(zoom)), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE, frame.plot = FALSE)
  axis(1, at = seq(from = start(zoom), to = end(zoom), by = (end(zoom) - start(zoom)) / 5), las = 1)

  mtext(ROI, outer = TRUE, cex = 1, line = -0.5)

  if (pdf | png) {
    dev.off()
  }
}

# plot4Cs_fig2B_reps(a=6, C17TetO, ChIPinfo = ChIPdf, domain_HMM_GR, zoomWidth = 4e6, yMax4C=1000, BWnorm='meanBG', outF = "/storage/shared/TACL/PK/manuscript/revision/fig2_revision/2B/pdf_reps/", scaleTriangle = 1e-08, pdf = FALSE, png=TRUE)


numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

foreach(a = seq_along(C17TetO), .packages = c("GenomicRanges", "graphics", "grDevices", "grid")) %dopar% {
   plot4Cs_fig2B_reps(a, C17TetO, ChIPinfo = ChIPdf, domain_HMM_GR, zoomWidth = 4e6, yMax4C=1000, BWnorm='meanBG', outF = paste0(outF, "/2B/pdf_reps_lwd1/"), scaleTriangle = 1e-08, pdf = TRUE, BWlwd =1)
}

# Stop the cluster
stopCluster(cl)




# 2C Tornado FLAG peaks within and outside domains --------------------------------

#dFLAG peaks.
#421 >35 & TACL HMM domain peaks


FLAG_filtered$dFLAG <- countOverlaps(FLAG_filtered,dFLAG)

expNames <- c(
    "C17_TACL-ON_FLAG_rep1_rep2", "C17_TACL-OFF_FLAG_rep2", "C17_mCherry_FLAG_rep1",
    "C17_TACL-ON_MAU2_rep1_rep2", "C17_TACL-OFF_MAU2_rep2", "C17_mCherry_MAU2_rep1_rep2",
    "C17_TACL-ON_NIPBL_rep2", "C17_TACL-OFF_NIPBL_rep2", "C17_mCherry_NIPBL_rep1_rep2",
    "C17_TACL-ON_RAD21_rep1_rep3", "C17_TACL-OFF_RAD21_rep3", "C17_mCherry_RAD21_rep1_rep3",
    "C17_TACL-ON_SMC1_rep2_rep3", "C17_TACL-OFF_SMC1_rep2_rep3", "C17_mCherry_SMC1_rep2_rep3",
    "C17_TACL-ON_STAG1_rep1", "C17_TACL-OFF_STAG1_rep1", "C17_mCherry_STAG1_rep1", "C17_WT_STAG1_rep1",
    "C17_TACL-ON_STAG2_rep1", "C17_TACL-OFF_STAG2_rep1", "C17_mCherry_STAG2_rep1",
    "C17_TACL-ON_CTCF_rep1", "C17_TACL-OFF_CTCF_rep1", "C17_mCherry_CTCF_rep1",
    "C17_TACL-ON_PDS5A_rep1", "C17_mCherry_PDS5A_rep1",
    "C17_TACL-ON_WAPL_rep1", "C17_mCherry_WAPL_rep1",
    "C17_TACL-ON_H3K27ac_rep1_rep2", "C17_TACL-OFF_H3K27ac_rep1_rep2", "C17_mCherry_H3K27ac_rep1_rep2"
)


ChIPdf_OI <- ChIPdf[ChIPdf$name %in% expNames, ]

toRnadoCov <- GetTornadoCov_V4(
  region = FLAG_filtered,
  ChIPdf = ChIPdf_OI,
  orderCov = "C17_TACL-ON_FLAG_rep1_rep2",
  orderby = "C17_TACL-ON_FLAG_rep1_rep2",
  zoom = 5000, binWidth = 10, Matrix = FALSE, makeList = FALSE, bincoverage = TRUE, prefix = prefix, bigwig = "pval"
)

saveRDS(toRnadoCov, file = file.path(covF, "Fig2C_dFLAG_pval_250224.rds"))

# groups
toRnadoCov <- readRDS(file.path(covF, "Fig2C_dFLAG_pval_250224.rds"))
toRnadoCov$peaks$dFLAG <- countOverlaps(toRnadoCov$peaks, dFLAG)
#table(toRnadoCov$peaks$dFLAG > 0, toRnadoCov$peaks$HMMdomain > 0) #421



toRnadoCov$peaks$group <- NA
toRnadoCov$peaks[toRnadoCov$peaks$distance > 3e6 & toRnadoCov$peaks$HMMdomain == 0]$group <- "GW" # 4795
toRnadoCov$peaks[toRnadoCov$peaks$HMMdomain > 0 & toRnadoCov$peaks$dFLAG > 0 ]$group <- "HMM dFLAG" #421
GROUPSdf <- as.data.frame(table(toRnadoCov$peaks$group))
plotGroups <- c("HMM dFLAG", "GW")

# groups that will be plotted
GROUPSdf[GROUPSdf$Var1 %in% plotGroups, ]



# order bwOI per protein

bwOI_df <-toRnadoCov$ChIPdf[c("name","condition", "prot", "rep")]

# Function to select the preferred replicate for each condition and protein
select_replicate <- function(df) {
  # Check if combined replicates are available
  combined_indices <- grep("rep1_rep2|rep2_rep3|rep1_rep3", df$rep)
  if (length(combined_indices) > 0) {
    return(df[combined_indices, , drop = FALSE]) # Return combined replicates
  } else {
    return(df[1, , drop = FALSE]) # Return the first available replicate if no combined
  }
}

# Split the data frame by condition and protein and apply function
split_df <- split(bwOI_df, list(bwOI_df$condition, bwOI_df$prot))
selected_entries <- lapply(split_df, select_replicate)

# Combine the results back into a single data frame
result_df <- do.call(rbind, selected_entries)
result_df <- result_df[!is.na(result_df$name),]

# Reorder the data frame by protein and condition
protein_order <- c("FLAG", "MAU2", "NIPBL", "SMC1", "RAD21", "STAG1", "STAG2", "CTCF", "PDS5A", "WAPL", "H3K27ac")
condition_order <- c("TACL-ON","TACL-OFF", "mCherry","WT")
result_df$prot <- factor(result_df$prot, levels = protein_order)
result_df$condition <- factor(result_df$condition, levels = condition_order)
bwOI_df_ordered <- result_df[order(result_df$prot, result_df$condition), ]
plotOrder <- bwOI_df_ordered$name


makePlotToRnado_V8(
  Plotregion_cov = toRnadoCov,
  plotOrder = plotOrder,
  orderFactor = toRnadoCov$infoDF$orderby[1],
  plotGroups = c("HMM dFLAG", "GW"),
  norm = "meanBG",
  normGroup = plotGroups,
  cutoffGroups = plotGroups,
  PlotMaxasCutoff = TRUE,
  SetCutOff = NULL,
  cutoffRatio = 0.8,
  cutoffLine = TRUE,
  subSampleGW = TRUE,
  subSampleGroups = "GW",
  # subsampleN=100,
  PeaksignalValue = 35,
  equalPlotMax = FALSE,
  setPlotMax = NULL,
  FlipRvStrand = FALSE,
  BigMatrix = TRUE,
  middleLine = FALSE,
  plotWidth = 2500,
  protColors = protCols,
  PNG = FALSE,
  PDF = TRUE,
  outF = file.path(outF, "2C/"),
  plotName = "2C_reps_GWFLAG_HMMdFLAG_pval_sigVal35_tornado_250224_ratio08",
  ColWidth = 4
)



# Fig 2D/E/F :ASP --------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


#ASP plots for each protein and condition ----------------------------------------------------------------------------------------------------


expNames <- c(
    "C17_TACL-ON_FLAG_rep1_rep2", "C17_TACL-OFF_FLAG_rep2",
    "C17_TACL-ON_MAU2_rep1_rep2", "C17_TACL-OFF_MAU2_rep2", "C17_mCherry_MAU2_rep1_rep2",
    "C17_TACL-ON_NIPBL_rep2", "C17_TACL-OFF_NIPBL_rep2", "C17_mCherry_NIPBL_rep1_rep2",
    "C17_TACL-ON_RAD21_rep1_rep3", "C17_TACL-OFF_RAD21_rep3", "C17_mCherry_RAD21_rep1_rep3",
    "C17_TACL-ON_SMC1_rep2_rep3", "C17_TACL-OFF_SMC1_rep2_rep3", "C17_mCherry_SMC1_rep2_rep3",
    "C17_TACL-ON_STAG1_rep1", "C17_TACL-OFF_STAG1_rep1", "C17_mCherry_STAG1_rep1",
    "C17_TACL-ON_STAG2_rep1", "C17_TACL-OFF_STAG2_rep1", "C17_mCherry_STAG2_rep1",
    "C17_TACL-ON_CTCF_rep1", "C17_TACL-OFF_CTCF_rep1", "C17_mCherry_CTCF_rep1",
    "C17_TACL-ON_PDS5A_rep1", "C17_mCherry_PDS5A_rep1",
    "C17_TACL-ON_WAPL_rep1", "C17_mCherry_WAPL_rep1",
    "C17_TACL-ON_H3K27ac_rep1_rep2", "C17_TACL-OFF_H3K27ac_rep1_rep2", "C17_mCherry_H3K27ac_rep1_rep2"
)
ChIP_factors <- c("FLAG", "MAU2", "NIPBL", "SMC1", "RAD21", "STAG1", "STAG2", "CTCF", "H3K27ac", "PDS5A", "WAPL")
ChIPdf_OI <- ChIPdf[ChIPdf$name %in% expNames, ]

toRnadoCov <- readRDS(file.path(covF, "Fig2C_dFLAG_pval_250224.rds"))
toRnadoCov$peaks$dFLAG <- countOverlaps(toRnadoCov$peaks, dFLAG)

toRnadoCov$peaks$group <- NA
toRnadoCov$peaks[toRnadoCov$peaks$distance > 3e6 & toRnadoCov$peaks$HMMdomain == 0]$group <- "GW" # 4795
toRnadoCov$peaks[toRnadoCov$peaks$HMMdomain > 0 & toRnadoCov$peaks$dFLAG > 0]$group <- "HMM dFLAG" # 421

plotGroups <- c("HMM dFLAG", "GW")

plotZoom <- 5000
plotWidth <- 1000
binWidth <- 10

# V151-V350
avgdf <- avgPlot_V2(toRnadoCov,
    plotOrder = expNames,
    #plotOrder = plotOrder[grepl(pattern = paste(ChIP_factors, collapse = "|"), plotOrder)],
    plotGroups = plotGroups,
    orderFactor = NULL,
    norm = "meanBG",
    #normGroup = plotGroups,
    PeaksignalValue = 0,
    plotWidth = plotWidth,
    FlipRvStrand = FALSE,
    middleLine = TRUE
)



middleBins <- c(plotZoom / binWidth / 2, (plotZoom / binWidth / 2) + 1)

avgdf_scaled <- avgdf

# Filter the data for the HMM dFLAG group
filtered_df <- avgdf_scaled %>% filter(group == "HMM dFLAG")

# Reshape the data to long format
filtered_df_long <- filtered_df %>%
  pivot_longer(cols = starts_with("V"), names_to = "Position", values_to = "Value") %>%
  mutate(Position = as.numeric(gsub("V", "", Position)))

# Calculate new positions based on bin size
filtered_df_long <- filtered_df_long %>%
  mutate(New_Position = (Position - middleBins[1]) * binWidth)

# Extract protein name from expName
filtered_df_long <- filtered_df_long %>%
  mutate(Protein = case_when(
    grepl("FLAG", expName) ~ "FLAG",
    grepl("MAU2", expName) ~ "MAU2",
    grepl("NIPBL", expName) ~ "NIPBL",
    grepl("SMC1", expName) ~ "SMC1",
    grepl("RAD21", expName) ~ "RAD21",
    grepl("STAG1", expName) ~ "STAG1",
    grepl("STAG2", expName) ~ "STAG2",
    grepl("CTCF", expName) ~ "CTCF",
    grepl("H3K27ac", expName) ~ "H3K27ac",
    grepl("PDS5A", expName) ~ "PDS5A",
    grepl("WAPL", expName) ~ "WAPL"
  ))


# Create a color mapping for TACL-ON, TACL-OFF, and mCherry
filtered_df_long <- filtered_df_long %>%
  #mutate(Condition = ifelse(grepl("ON", expName), "TACL-ON", "TACL-OFF"))
  mutate(Condition = case_when(
    grepl("ON", expName) ~ "TACL-ON",
    grepl("OFF", expName) ~ "TACL-OFF",
    grepl("mCherry", expName) ~ "mCherry"
  ))


y_ranges <- filtered_df_long %>%
  group_by(Protein) %>%
  summarize(
    y_min = 0,
    y_max = ceiling(max(Value, na.rm = TRUE) * 2) / 2  # round up to nearest 0.5
  )
)



# Create a line plot for each protein
proteins <- unique(filtered_df_long$Protein)

plot_list <- list()


for (protein in proteins) {
  plot_data <- filtered_df_long %>% filter(Protein == protein)
  y_range <- y_ranges %>% filter(Protein == protein)

  p <- ggplot(plot_data, aes(x = New_Position, y = Value, group = expName, color = Condition)) +
  
    geom_line(linewidth = 1.2) +
    labs(
      title = protein,
      x = "Position (bp)",
      y = "Normalized Signal",
      color = "Condition"
    ) +
    scale_y_continuous(
      limits = c(0, y_range$y_max),
      expand = c(0, 0),
      breaks = seq(0, y_range$y_max, 0.5) # Set breaks at every integer
    ) +
    scale_color_manual(values = c("TACL-ON" = "#1f78b4", "TACL-OFF" = "#e31a1c", "mCherry" = "#33a02c")) +
    guides(color = guide_legend(override.aes = list(group = NULL))) + 
    theme_minimal(base_size = 18) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, margin = margin(r = 10)),
      axis.text = element_text(size = 14),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14)
    )

  # Store each plot in a list using protein names as keys
  plot_list[[protein]] <- p
}

# Combine all plots using patchwork with a shared legend and title
combined_plot <- wrap_plots(plot_list) +
  plot_layout(guides = 'collect') &
  plot_annotation(title = "Normalize Signal at Local dFLAG Peaks") &
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.position = 'bottom'
  )

ggsave(file.path(outF, "2D/All_ASP_combined_scaledGW_250327.pdf"), combined_plot, width = 15, height = 8)




# Fig 2H : Tornado CTCF+/dFLAG peaks inside TACL domains --------------------------------



table(dFLAG$domain > 0) # 454
dFLAG$filteredFlag <- countOverlaps(dFLAG, FLAG_filtered)
table(dFLAG$filteredFlag > 0) # 460
table(dFLAG$filteredFlag > 0, dFLAG$domain > 0) # 421  --421 dFLAG filter peaks within the HMM domain

local_filtered_dFLAG <- dFLAG[dFLAG$filteredFlag > 0 & dFLAG$domain > 0] # 421 peaks
local_filtered_dFLAG$type <- "dFLAG"

#FLAG_filtered$dFLAG <- countOverlaps(FLAG_filtered,dFLAG)
#table(FLAG_filtered$dFLAG > 0, FLAG_filtered$signalValue>35) 


CTCF_peaks <- peak_CTCF_C17_TMAU2[peak_CTCF_C17_TMAU2$signalValue > 35] # 27941
CTCF_peaks$type <- "CTCF"
CTCF_peaks$domain <- CTCF_peaks$HMMdomain

CTCF_peaks$dFLAG <- countOverlaps(CTCF_peaks, dFLAG) # 421
CTCF_peaks$local_filtered_dFLAG <- countOverlaps(CTCF_peaks, local_filtered_dFLAG) # 421

table(CTCF_peaks$dFLAG > 0, CTCF_peaks$HMMdomain > 0) # 246 OL with dFLAG, 366 not
table(CTCF_peaks$local_filtered_dFLAG > 0, CTCF_peaks$HMMdomain > 0) #239 OL with filtered dFLAG, 373 not

#combine the peaks of interest (dFLAG and CTCF) to make the tornado plot
CTCF_dFLAG_peaks <- c(local_filtered_dFLAG[,c("distance","domain","type")], 
                      CTCF_peaks[,c("distance","domain","type")]
                      )

table(CTCF_dFLAG_peaks$type) # 421 dFLAG, 27941 CTCF






expNames <- c(
    "C17_TACL-ON_FLAG_rep1_rep2", "C17_TACL-OFF_FLAG_rep2", "C17_mCherry_FLAG_rep1",
    "C17_TACL-ON_RAD21_rep1_rep3", "C17_TACL-OFF_RAD21_rep3", "C17_mCherry_RAD21_rep1_rep3",
    "C17_TACL-ON_SMC1_rep2_rep3", "C17_TACL-OFF_SMC1_rep2_rep3", "C17_mCherry_SMC1_rep2_rep3",
    "C17_TACL-ON_STAG1_rep1", "C17_TACL-OFF_STAG1_rep1", "C17_mCherry_STAG1_rep1",
    "C17_TACL-ON_STAG2_rep1", "C17_TACL-OFF_STAG2_rep1", "C17_mCherry_STAG2_rep1",
    "C17_TACL-ON_CTCF_rep1", "C17_TACL-OFF_CTCF_rep1", "C17_mCherry_CTCF_rep1",
    "C17_TACL-ON_PDS5A_rep1", "C17_mCherry_PDS5A_rep1",
    "C17_TACL-ON_WAPL_rep1", "C17_mCherry_WAPL_rep1"
)



ChIPdf_OI <- ChIPdf[ChIPdf$name %in% expNames, ]

toRnadoCov <- GetTornadoCov_V4(
  region = CTCF_dFLAG_peaks,
  ChIPdf = ChIPdf_OI,
  orderCov = "C17_TACL-ON_FLAG_rep1_rep2",
  orderby = "C17_TACL-ON_FLAG_rep1_rep2",
  zoom = 5000, binWidth = 10, Matrix = FALSE, makeList = FALSE, bincoverage = TRUE, prefix = prefix, bigwig = "pval"
)

saveRDS(toRnadoCov, file = file.path(covF, "Fig2C_CTCF_dFLAG_pval_250330.rds"))




# groups
toRnadoCov <- readRDS(file.path(covF, "Fig2C_CTCF_dFLAG_pval_250330.rds"))


toRnadoCov$peaks$group <- NA
toRnadoCov$peaks[which(toRnadoCov$peaks$distance > 3e6 & toRnadoCov$peaks$domain == 0 & toRnadoCov$peaks$type == "CTCF")]$group <- "GW_CTCF" # 4795
toRnadoCov$peaks[which(toRnadoCov$peaks$domain > 0 & toRnadoCov$peaks$type == "dFLAG") ]$group <- "HMM_dFLAG" #421
GROUPSdf <- as.data.frame(table(toRnadoCov$peaks$group))
plotGroups <- c("HMM_dFLAG", "GW_CTCF")

# groups that will be plotted
GROUPSdf[GROUPSdf$Var1 %in% plotGroups, ]



# order bwOI per protein

bwOI_df <-toRnadoCov$ChIPdf[c("name","condition", "prot", "rep")]

# Function to select the preferred replicate for each condition and protein
select_replicate <- function(df) {
  # Check if combined replicates are available
  combined_indices <- grep("rep1_rep2|rep2_rep3|rep1_rep3", df$rep)
  if (length(combined_indices) > 0) {
    return(df[combined_indices, , drop = FALSE]) # Return combined replicates
  } else {
    return(df[1, , drop = FALSE]) # Return the first available replicate if no combined
  }
}

# Split the data frame by condition and protein and apply function
split_df <- split(bwOI_df, list(bwOI_df$condition, bwOI_df$prot))
selected_entries <- lapply(split_df, select_replicate)

# Combine the results back into a single data frame
result_df <- do.call(rbind, selected_entries)
result_df <- result_df[!is.na(result_df$name),]

# Reorder the data frame by protein and condition
protein_order <- c("FLAG", "SMC1", "RAD21", "STAG1", "STAG2", "CTCF", "PDS5A", "WAPL")
condition_order <- c("TACL-ON","TACL-OFF", "mCherry")
result_df$prot <- factor(result_df$prot, levels = protein_order)
result_df$condition <- factor(result_df$condition, levels = condition_order)
bwOI_df_ordered <- result_df[order(result_df$prot, result_df$condition), ]
plotOrder <- bwOI_df_ordered$name


makePlotToRnado_V8(
  Plotregion_cov = toRnadoCov,
  plotOrder = plotOrder,
  orderFactor = toRnadoCov$infoDF$orderby[1],
  plotGroups = plotGroups,
  norm = "meanBG",
  normGroup = plotGroups,
  cutoffGroups = plotGroups,
  PlotMaxasCutoff = TRUE,
  SetCutOff = NULL,
  cutoffRatio = 0.8,
  cutoffLine = TRUE,
  subSampleGW = TRUE,
  subSampleGroups = "GW_CTCF",
  # subsampleN=100,
  PeaksignalValue = 0,  #alerady prefiltered for >35
  equalPlotMax = FALSE,
  setPlotMax = NULL,
  FlipRvStrand = FALSE,
  BigMatrix = TRUE,
  middleLine = FALSE,
  plotWidth = 2500,
  protColors = protCols,
  PNG = FALSE,
  PDF = TRUE,
  outF = file.path(outF, "2C/"),
  plotName = "2C_GWCTCF_HMMdFLAG_pval_sigVal35_tornado_250330_ratio08",
  ColWidth = 4
)


