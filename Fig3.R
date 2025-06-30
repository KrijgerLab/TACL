

#Figure 3 ------------------------------------

# functions ----------------------------------------------------------------

# Add all functions to source file


prefix <- ""
source(file.path(prefix, "./TACL_functions_NG2025.r"))
source(file.path(prefix, "./TACL_features_NG2025.r"))


outF <- "./out/Fig3"
if (!dir.exists(outF)) {
  dir.create(outF, recursive = TRUE)
}


#Fig 3A ----------------------------------------


#NIPBL & SMC1 peaks in RAD21 AID clone 
 NIPBL_RAD21AID_opt<-import(paste0(prefix,'/TACL/data/peaks/VER6935_TACL_RH_EHap1_C17_TetR-hMAU2v2_Rad21-AID-GB_NIPBL_spikein_optimal_peak.regionPeak.bb'))
 NIPBL_RAD21AID_opt<-NIPBL_RAD21AID_opt[NIPBL_RAD21AID_opt$signalValue>= 6]

 NIPBL_RAD21AID_IAA_opt<-import(paste0(prefix,'/TACL/data/peaks/VER6935_TACL_RH_EHap1_C17_TetR-hMAU2v2_Rad21-AID-GB_IAA2h_NIPBL_spikein_optimal_peak.regionPeak.bb'))
 NIPBL_RAD21AID_IAA_opt<-NIPBL_RAD21AID_IAA_opt[NIPBL_RAD21AID_IAA_opt$signalValue>= 4]

 NIPBL_RAD21AID_combined<-c(NIPBL_RAD21AID_opt,NIPBL_RAD21AID_IAA_opt)
 NIPBL_RAD21AID_combined<-filterOverlapPeaks(peaks=NIPBL_RAD21AID_combined)


 SMC1_RAD21AID_opt<-import(paste0(prefix,'/TACL/data/peaks/VER6935_TACL_RH_EHap1_C17_TetR-hMAU2v2_Rad21-AID-GB_SMC1_spikein_optimal_peak.regionPeak.bb'))
 SMC1_RAD21AID_opt<-SMC1_RAD21AID_opt[SMC1_RAD21AID_opt$signalValue>= 15]

 SMC1_RAD21AID_IAA_opt<-import(paste0(prefix,'/TACL/data/peaks/VER6935_TACL_RH_EHap1_C17_TetR-hMAU2v2_Rad21-AID-GB_IAA2h_SMC1_spikein_optimal_peak.regionPeak.bb'))
 SMC1_RAD21AID_IAA_opt<-SMC1_RAD21AID_IAA_opt[SMC1_RAD21AID_IAA_opt$signalValue>= 10]

 SMC1_RAD21AID_combined<-c(SMC1_RAD21AID_opt,SMC1_RAD21AID_IAA_opt)
 SMC1_RAD21AID_combined<-filterOverlapPeaks(peaks=SMC1_RAD21AID_combined)

 SMC1_NIPBL_combinedPeaks<-c(NIPBL_RAD21AID_combined,SMC1_RAD21AID_combined) #22204
 SMC1_NIPBL_combinedPeaks<-GenomicRanges::reduce(SMC1_NIPBL_combinedPeaks) #21292 (912 OL)

 SMC1_NIPBL_combinedPeaks<-addTetOdistance(SMC1_NIPBL_combinedPeaks,C17TetO)
 SMC1_NIPBL_combinedPeaks$signalValue<-1 #placeholder as after reduction we no longer have this.

 SMC1_NIPBL_combinedPeaks$NIPBL<-countOverlaps(query = SMC1_NIPBL_combinedPeaks, NIPBL_RAD21AID_combined) #1830
 SMC1_NIPBL_combinedPeaks$SMC1<-countOverlaps(query = SMC1_NIPBL_combinedPeaks, SMC1_RAD21AID_combined) #20366 
 




 #coverage
expNames_A <- ChIPdf[
                ChIPdf$cell %in% c("C17_RAD21AID_SC1","C17_RAD21AID_D8", "C17_RAD21AID_SC6") & 
                ChIPdf$prot %in% c("NIPBL", "MAU2") &
                ChIPdf$condition %in% c("TACL-ON", "TACL-ON_IAA", "mCherry", "mCherry_IAA")
                ,"name"]

expNames_B <- ChIPdf[
                ChIPdf$cell %in% c("C17", "HAP1") & 
                ChIPdf$prot %in% c("H3K4me1", "CTCF") &
                ChIPdf$condition %in% c("TACL-ON", "mCherry", "WT")
                ,"name"]
#"C17_mCherry_CTCF_rep1"     "HAP1_WT_H3K4me1_rep1_rep2"

expNames_C <- ChIPdf[
                ChIPdf$cell %in% c("C17", "HAP1") & 
                ChIPdf$prot %in% c("H3K27ac", "H3K4me3") &
                ChIPdf$condition %in% c("TACL-ON","mCherry") &
                ChIPdf$rep %in% c ("rep1_rep2")
                ,"name"]
#"C17_mCherry_H3K27ac_rep1_rep2" "C17_mCherry_H3K4me3_rep1_rep2"


expNames_D <- ChIPdf[
                ChIPdf$cell %in% c("C17_RAD21AID_SC1","C17_RAD21AID_D8", "C17_RAD21AID_SC6") & 
                ChIPdf$prot %in% c("RAD21", "SMC1", "SMC3") &
                ChIPdf$condition %in% c("TACL-ON", "TACL-ON_IAA", "mCherry", "mCherry_IAA")
                ,"name"]


expNames <- c(expNames_A, expNames_B, expNames_C, expNames_D)
ChIPdf_OI <- ChIPdf[ChIPdf$name %in% expNames, ]


toRnadoCov <- GetTornadoCov_V4(
                region = SMC1_NIPBL_combinedPeaks,
                ChIPdf = ChIPdf_OI,
                orderCov = expNames[1:8],
                orderby = "C17_RAD21AID_SC1_TACL-ON_IAA_NIPBL_rep1",
                zoom = 5000, binWidth = 10, Matrix = FALSE, makeList = FALSE, bincoverage = TRUE, prefix = prefix, bigwig='pval'
            )

saveRDS(toRnadoCov, file = file.path(covF, "3A_SMC1NIPBL_RAD21SMC1SMC3_pval_240702.rds"))




toRnadoCov$peaks$domain <- countOverlaps(toRnadoCov$peaks, domain_HMM_GR)
toRnadoCov$peaks <- addTetOdistance(toRnadoCov$peaks, C17TetO)
toRnadoCov$peaks$CTCF <- countOverlaps(toRnadoCov$peaks, CTCF_reduced)
toRnadoCov$peaks$FLAG <- countOverlaps(toRnadoCov$peaks, FLAG_filtered)
toRnadoCov$peaks$H3K27ac <- countOverlaps(toRnadoCov$peaks, NKI_H3K27ac_peaks)
toRnadoCov$peaks$H3K4me1 <- countOverlaps(toRnadoCov$peaks, H3K4me1_opt_reduced)
toRnadoCov$peaks$H3K4me3 <- countOverlaps(toRnadoCov$peaks, NKI_H3K4me3_peaks)


#GW Enh; Prom; CTCF
toRnadoCov$peaks$group <- NA
toRnadoCov$peaks[which(toRnadoCov$peaks$distance > 3e6 & toRnadoCov$peaks$domain == 0 & toRnadoCov$peaks$H3K4me3 ==0 & toRnadoCov$peaks$H3K4me1>0) ]$group <-"GW enhancer"
toRnadoCov$peaks[which(toRnadoCov$peaks$distance > 3e6 & toRnadoCov$peaks$domain == 0 & toRnadoCov$peaks$H3K4me3 ==0 & toRnadoCov$peaks$H3K27ac>0) ]$group <-"GW enhancer"

toRnadoCov$peaks[which(toRnadoCov$peaks$distance > 3e6 & toRnadoCov$peaks$domain == 0 & toRnadoCov$peaks$H3K4me3 > 0) ]$group <-"GW promoter"

toRnadoCov$peaks[which(toRnadoCov$peaks$distance > 3e6 & toRnadoCov$peaks$domain == 0 & toRnadoCov$peaks$H3K4me3 ==0 & toRnadoCov$peaks$H3K4me1==0 & toRnadoCov$peaks$H3K27ac==0 & toRnadoCov$peaks$CTCF==0) ]$group <-"GW prom-;enh-; CTCF-"


toRnadoCov$peaks[which(toRnadoCov$peaks$distance > 3e6 & toRnadoCov$peaks$domain == 0 & toRnadoCov$peaks$CTCF > 0) ]$group <-"GW CTCF"



GROUPSdf <- as.data.frame(table(toRnadoCov$peaks$group))

#sampleGroups <- GROUPSdf[GROUPSdf$Freq > 1000, "Var1"]


plotGroups <- GROUPSdf$Var1

# groups that will be plotted
GROUPSdf[GROUPSdf$Var1 %in% plotGroups, ]

# which groups will not be plotted
GROUPSdf[!GROUPSdf$Var1 %in% plotGroups, ]

# order bwOI per protein


bwOI_df <- data.frame(name = expNames, prot = ChIPdf[match(expNames, ChIPdf$name), "prot"], condition = ChIPdf[match(expNames, ChIPdf$name), "condition"])



prot_order <- c("SMC1", "RAD21","SMC3","MAU2","NIPBL", "CTCF", "H3K4me1", "H3K27ac", "H3K4me3")
bwOI_df$prot <- factor(bwOI_df$prot, levels = prot_order)
condition_order <- c("TACL-ON", "TACL-ON_IAA","mCherry", "mCherry_IAA","WT")
bwOI_df$condition <- factor(bwOI_df$condition, levels = condition_order)
bwOI_df_ordered <- bwOI_df[order(bwOI_df$prot, bwOI_df$condition), ]
plotOrder <- bwOI_df_ordered$name
any(plotOrder %in% toRnadoCov$ChIPdf$name)

plotOrder <- c(
    bwOI_df_ordered$name
)


# > plotOrder
#  [1] "C17_RAD21AID_SC1_TACL-ON_SMC1_rep1"          "C17_RAD21AID_SC1_TACL-ON_SMC1_rep2"          "C17_RAD21AID_SC1_TACL-ON_SMC1_rep1_rep2"
#  [4] "C17_RAD21AID_SC1_TACL-ON_IAA_SMC1_rep1"      "C17_RAD21AID_SC1_TACL-ON_IAA_SMC1_rep2"      "C17_RAD21AID_SC1_TACL-ON_IAA_SMC1_rep1_rep2"
#  [7] "C17_RAD21AID_SC6_mCherry_SMC1_rep1"          "C17_RAD21AID_SC6_mCherry_SMC1_rep2"          "C17_RAD21AID_SC6_mCherry_SMC1_rep1_rep2"
# [10] "C17_RAD21AID_SC6_mCherry_IAA_SMC1_rep1"      "C17_RAD21AID_SC1_TACL-ON_RAD21_rep1"         "C17_RAD21AID_SC1_TACL-ON_IAA_RAD21_rep1"
# [13] "C17_RAD21AID_SC6_mCherry_RAD21_rep1"         "C17_RAD21AID_SC6_mCherry_IAA_RAD21_rep1"     "C17_RAD21AID_SC1_TACL-ON_SMC3_rep1"
# [16] "C17_RAD21AID_SC1_TACL-ON_IAA_SMC3_rep1"      "C17_RAD21AID_SC6_mCherry_SMC3_rep1"          "C17_RAD21AID_SC1_TACL-ON_MAU2_rep1"
# [19] "C17_RAD21AID_SC1_TACL-ON_IAA_MAU2_rep1"      "C17_RAD21AID_D8_mCherry_MAU2_rep1"           "C17_RAD21AID_D8_mCherry_IAA_MAU2_rep1"
# [22] "C17_RAD21AID_SC1_TACL-ON_NIPBL_rep1"         "C17_RAD21AID_SC1_TACL-ON_IAA_NIPBL_rep1"     "C17_RAD21AID_SC6_mCherry_NIPBL_rep1"
# [25] "C17_RAD21AID_SC6_mCherry_IAA_NIPBL_rep1"     "C17_mCherry_CTCF_rep1"                       "HAP1_WT_H3K4me1_rep1_rep2"
# [28] "C17_mCherry_H3K27ac_rep1_rep2"               "C17_mCherry_H3K4me3_rep1_rep2"


plotOrder_TACL <- grep("TACL-ON", plotOrder, value = TRUE)
plotOrder_TACL_REPS <- c("C17_RAD21AID_SC1_TACL-ON_SMC1_rep1_rep2", "C17_RAD21AID_SC1_TACL-ON_IAA_SMC1_rep1_rep2", plotOrder_TACL[7:17], "HAP1_WT_H3K4me1_rep1_rep2")

plotOrder_mCherry <- grep("mCherry", plotOrder, value = TRUE)
plotOrder_mCherry_REPS <- c("C17_RAD21AID_SC6_mCherry_SMC1_rep1_rep2", plotOrder_mCherry[c(4:6, 8:14)], "HAP1_WT_H3K4me1_rep1_rep2")


# order bwOI per protein
#any(plotOrder %in% toRnadoCov$ChIPdf$name)



makePlotToRnado_V8(
    Plotregion_cov = toRnadoCov,
    plotOrder = plotOrder_TACL_REPS,
    orderFactor = "C17_RAD21AID_SC1_TACL-ON_IAA_NIPBL_rep1",
    plotGroups = plotGroups,
    norm = "outerBG",
    #normGroup = plotGroups[grepl("GW", plotGroups)],
    PlotMaxasCutoff = TRUE,
    SetCutOff = NULL,
    cutoffRatio = 0.8,
    cutoffLine = TRUE,
    subSampleGroups = GROUPSdf[GROUPSdf$Freq > 1000, "Var1"],
    subSampleGW = TRUE,
    #subsampleN=1000,
    PeaksignalValue = 0,
    equalPlotMax = FALSE,
    setPlotMax = NULL,
    FlipRvStrand = FALSE,
    BigMatrix = TRUE,
    middleLine = FALSE,
    plotWidth = 2500,
    protColors = protCols,
    useRaster = TRUE,
    PNG = FALSE,
    PDF = TRUE,
    outF = outF,
    plotName = "3A_TACL-ON_SMC1NIPBL_CTCFPromEnh_pval_240701_ratio08",
    ColWidth = 4
)



#mCherry plot
makePlotToRnado_V8(
  Plotregion_cov = toRnadoCov,
  plotOrder = plotOrder_mCherry_REPS,
  orderFactor = "C17_RAD21AID_SC1_TACL-ON_IAA_NIPBL_rep1",
  plotGroups = plotGroups,
  norm = "outerBG",
  # normGroup = plotGroups[grepl("GW", plotGroups)],
  PlotMaxasCutoff = TRUE,
  SetCutOff = NULL,
  cutoffRatio = 0.8,
  cutoffLine = TRUE,
  subSampleGroups = GROUPSdf[GROUPSdf$Freq > 1000, "Var1"],
  subSampleGW = TRUE,
  # subsampleN=1000,
  PeaksignalValue = 0,
  equalPlotMax = FALSE,
  setPlotMax = NULL,
  FlipRvStrand = FALSE,
  BigMatrix = TRUE,
  middleLine = FALSE,
  plotWidth = 2500,
  protColors = protCols,
  useRaster = TRUE,
  PNG = FALSE,
  PDF = TRUE,
  outF = outF,
  plotName = "3A_mCherry_SMC1NIPBL_CTCFPromEnh_pval_240701_ratio08",
  ColWidth = 4
)



# 3D ------------------------------------------------------------
#3D_flag_pval_NIPBL-V5_FLAG_filtered_by_C17_V5-MAU2_TACL-ON_V5_rep1_0.8_FlipFALSE_outerBins2024-07-02.pdf



FLAG_filtered$CTCF <- countOverlaps(FLAG_filtered, CTCF_reduced)
FLAG_filtered$dFLAG <- countOverlaps(FLAG_filtered, dFLAG)

expNames_A <- ChIPdf[
  ChIPdf$cell %in% c("C17_V5-MAU2") &
    ChIPdf$prot %in% c("V5", "NIPBL") &
    ChIPdf$condition %in% c("mCherry", "TACL-ON")
  # & ChIPdf$rep %in% paste0("rep", 1:10)
  , "name"
]
expNames_B <- ChIPdf[
  ChIPdf$cell %in% c("C17") &
    ChIPdf$prot %in% c("FLAG", "CTCF", "NIPBL") &
    ChIPdf$condition %in% c("TACL-ON", "mCherry")
  # & ChIPdf$rep %in% paste0("rep", 1:10)
  , "name"
]


expNames_B <- c("C17_TACL-ON_CTCF_rep1", "C17_TACL-ON_NIPBL_rep2", "C17_TACL-ON_FLAG_rep1_rep2", "C17_mCherry_CTCF_rep1", "C17_mCherry_FLAG_rep1", "C17_mCherry_NIPBL_rep1_rep2")

expNames <- c(expNames_A, expNames_B)
ChIPdf_OI <- ChIPdf[ChIPdf$name %in% expNames, ]



toRnadoCov <- GetTornadoCov_V4(
  region = FLAG_filtered,
  ChIPdf = ChIPdf_OI,
  orderCov = expNames,
  orderby = "C17_V5-MAU2_TACL-ON_V5_rep1",
  zoom = 5000, binWidth = 10, Matrix = FALSE, makeList = FALSE, bincoverage = TRUE, prefix = prefix, bigwig = "pval"
)

saveRDS(toRnadoCov, file = file.path(covF, "3D_FLAGfilt_pval_240702.rds"))



# Add groups
toRnadoCov$peaks$group <- NA
toRnadoCov$peaks[which(toRnadoCov$peaks$distance > 3e6 & toRnadoCov$peaks$HMMdomain == 0)]$group <- "GW"
toRnadoCov$peaks[which(toRnadoCov$peaks$HMMdomain > 0 & toRnadoCov$peaks$dFLAG > 0 & toRnadoCov$peaks$CTCF > 0)]$group <- "Local dFLAG;CTCF+"
toRnadoCov$peaks[which(toRnadoCov$peaks$HMMdomain > 0 & toRnadoCov$peaks$dFLAG > 0 & toRnadoCov$peaks$CTCF == 0)]$group <- "Local dFLAG;CTCF-"


GROUPSdf <- as.data.frame(table(toRnadoCov$peaks$group))

# sampleGroups <- GROUPSdf[GROUPSdf$Freq > 1000, "Var1"]



plotGroups <- GROUPSdf$Var1

# groups that will be plotted
GROUPSdf[GROUPSdf$Var1 %in% plotGroups, ]

# which groups will not be plotted
GROUPSdf[!GROUPSdf$Var1 %in% plotGroups, ]

# order bwOI per protein


bwOI_df <- data.frame(name = expNames, prot = ChIPdf[match(expNames, ChIPdf$name), "prot"], condition = ChIPdf[match(expNames, ChIPdf$name), "condition"])




prot_order <- c("V5", "NIPBL", "FLAG", "CTCF")
bwOI_df$prot <- factor(bwOI_df$prot, levels = prot_order)
condition_order <- c("TACL-ON", "mCherry")
bwOI_df$condition <- factor(bwOI_df$condition, levels = condition_order)
bwOI_df_ordered <- bwOI_df[order(bwOI_df$condition, bwOI_df$prot), ]
plotOrder <- bwOI_df_ordered$name
any(plotOrder %in% toRnadoCov$ChIPdf$name)

plotOrder <- c(
  bwOI_df_ordered$name
)


# order bwOI per protein
any(plotOrder %in% toRnadoCov$ChIPdf$name)



makePlotToRnado_V8(
  Plotregion_cov = toRnadoCov,
  plotOrder = plotOrder,
  orderFactor = "C17_V5-MAU2_TACL-ON_V5_rep1",
  plotGroups = plotGroups,
  norm = "outerBins",
  normGroup = plotGroups,
  PlotMaxasCutoff = TRUE,
  SetCutOff = NULL,
  cutoffRatio = 0.8,
  cutoffLine = TRUE,
  subSampleGroups = GROUPSdf[GROUPSdf$Freq > 1000, "Var1"],
  subSampleGW = TRUE,
  # subsampleN=20000,
  PeaksignalValue = 35,
  equalPlotMax = FALSE,
  setPlotMax = NULL,
  FlipRvStrand = FALSE,
  BigMatrix = TRUE,
  middleLine = FALSE,
  plotWidth = 2500,
  protColors = protCols,
  useRaster = TRUE,
  PNG = FALSE,
  PDF = TRUE,
  outF = file.path(outF, "3D/"),
  plotName = "3D_flag_pval_NIPBL-V5",
  ColWidth = 4
)





#3E Violin ----------------------------------------


# Violin plot with log transformation ----------------------------------------------------------



FLAG_filtered$dFLAG <- countOverlaps(FLAG_filtered, dFLAG)
FLAG_filtered$CTCF <- countOverlaps(FLAG_filtered, CTCF_reduced)
FLAG_filtered$CTCF_on <- countOverlaps(FLAG_filtered, peak_CTCF_C17_TMAU2)
FLAG_filtered$CTCF_mch <- countOverlaps(FLAG_filtered, peak_CTCF_C17_Tmcherry)


expNames <- c("C17_TACL-ON_FLAG_rep1_rep2","C17_TACL-OFF_FLAG_rep2", "C17_TACL-ON_NIPBL_rep2", "C17_TACL-OFF_NIPBL_rep2", "C17_V5-MAU2_TACL-ON_V5_rep1", "C17_V5-MAU2_TACL-ON_NIPBL_rep1","C17_V5-MAU2_mCherry_NIPBL_rep1","C17_V5-MAU2_mCherry_V5_rep1")
ChIPdf_OI <- ChIPdf[ChIPdf$name %in% expNames, ]



region <- FLAG_filtered
for(a in 1:nrow(ChIPdf_OI)){
  bw <- paste0(prefix, ChIPdf_OI[a, "bigwig"])
  name <- ChIPdf_OI[a, "name"]
  message(name)
  region <- getCov(bw, region, name)
}


#scale the data

scaled_region <- region

for (i in seq_len(nrow(ChIPdf_OI))) {
  cov_col <- ChIPdf_OI$name[i]
  meanBG <- ChIPdf_OI$meanBG[i]

  if (cov_col %in% colnames(mcols(scaled_region))) {
    # Scale the coverage values by dividing by meanBG
    message(paste0("Scaling ", cov_col, " by ", meanBG))
    mcols(scaled_region)[[cov_col]] <- mcols(scaled_region)[[cov_col]] / meanBG
  }
}



scaled_region$group <- "other"
scaled_region$group[scaled_region$HMMdomain > 0] <- "local"
scaled_region$group[scaled_region$HMMdomain == 0 & scaled_region$distance>3e6] <- "GW"


scaled_region_filtered <- scaled_region[scaled_region$group %in% c("local", "GW"), ]
FLAG_CTCF_peaks <- scaled_region_filtered[scaled_region_filtered$CTCF > 0, ] #299
FLAG_NOCTCF_peaks <- scaled_region_filtered[scaled_region_filtered$CTCF == 0, ] #205


df <- as.data.frame(scaled_region_filtered)
df$group_condition <- with(df, ifelse(dFLAG == 0 & distance>3e6 & HMMdomain==0, "GW dFLAG == 0",
                                      ifelse(dFLAG > 0 & CTCF == 0 & HMMdomain>0, "Local dFLAG > 0 & CTCF == 0", 
                                             ifelse(dFLAG > 0 & CTCF > 0  & HMMdomain>0, "Local dFLAG > 0 & CTCF > 0", NA))))







# Filter and mutate data
df_filtered <- df[!is.na(df$group_condition), ]
df_filtered <- df_filtered[df_filtered$group_condition %in% c("GW dFLAG == 0", "Local dFLAG > 0 & CTCF > 0"), ]  #4766; 267

# Ensure all relevant columns are numeric and add pseudocount
numeric_columns <- c("C17_TACL.ON_FLAG_rep1_rep2", "C17_TACL.OFF_FLAG_rep2", 
                     "C17_V5.MAU2_TACL.ON_V5_rep1", "C17_V5.MAU2_TACL.ON_NIPBL_rep1", 
                     "C17_V5.MAU2_mCherry_NIPBL_rep1", "C17_V5.MAU2_mCherry_V5_rep1")
pseudocount <- 1e-6
df_filtered[numeric_columns] <- lapply(df_filtered[numeric_columns], function(x) as.numeric(x) + pseudocount)





# Create the violin plot with log transformation, point overlay, and customizations


library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# Calculate the fold change
df_filtered <- df_filtered %>%
  mutate(C17_TACL.ON_NIPBL_rep2_fc = C17_TACL.ON_NIPBL_rep2 / C17_TACL.OFF_NIPBL_rep2,
         C17_TACL.ON_FLAG_rep1_rep2_fc = C17_TACL.ON_FLAG_rep1_rep2 / C17_TACL.OFF_FLAG_rep2,
         C17_V5.MAU2_TACL.ON_NIPBL_rep1_fc = C17_V5.MAU2_TACL.ON_NIPBL_rep1 / C17_V5.MAU2_mCherry_NIPBL_rep1,
         C17_V5.MAU2_TACL.ON_V5_rep1_fc = C17_V5.MAU2_TACL.ON_V5_rep1 / C17_V5.MAU2_mCherry_V5_rep1)

# Apply log2 transformation to the fold change
df_filtered <- df_filtered %>%
  mutate(across(contains("fc"), log2))

# Compute y-axis limits for the fold change columns
fc_y_limits <- range(df_filtered %>% select(contains("fc")), na.rm = TRUE)

# Pivot the dataframe to long format for plotting
df_long <- df_filtered %>%
  pivot_longer(cols = contains("fc"), names_to = "variable", values_to = "value")


df_long$group_condition <- factor(df_long$group_condition, levels = c("GW dFLAG == 0", "Local dFLAG > 0 & CTCF > 0"))

custom_colors <- c(
  "GW dFLAG == 0" = "#E69F00",  
  "Local dFLAG > 0 & CTCF > 0" = "#56B4E9"
)


p_log2_fc <- ggplot(df_long, aes(x = variable, y = value, fill = group_condition)) +
  geom_violin(position = position_dodge(width = 0.9), alpha = 0.7) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Coverage at FLAG Peaks",
       x = "Condition",
       y = "Fold Change (log2)") +
  theme_minimal()

print(p_log2_fc)





ggsave(file.path(outF,"3D_violin.pdf"), p_log2_fc, width = 10, height = 7)




t_test_results <- df_long %>%
  group_by(variable) %>%
  summarize(
    p_value = t.test(value ~ group_condition)$p.value
  )

# Print the p-values
print(t_test_results)


# Perform Mann-Whitney U-test for each condition (variable) between the two groups
mann_whitney_results <- df_long %>%
  group_by(variable) %>%
  summarize(
    p_value = wilcox.test(value ~ group_condition)$p.value
  )

# Print the results
print(mann_whitney_results)





# 3F: Substraction Violin CTCF+/dFLAG  ----------------------------------------



outF <- "/storage/shared/TACL/PK/manuscript/revision/fig/2E/"
dir.create(outF, showWarnings = FALSE)




# CTCF peaks
peak_CTCF_C17_TMAU2$dFLAG <- countOverlaps(peak_CTCF_C17_TMAU2, dFLAG)
peak_CTCF_C17_TMAU2$domain <- countOverlaps(peak_CTCF_C17_TMAU2, domain_HMM_GR)
CTCF_dFLAG_HMM <- peak_CTCF_C17_TMAU2[peak_CTCF_C17_TMAU2$domain > 0 & peak_CTCF_C17_TMAU2$dFLAG > 0]

# Give the downstream TetO motifs a negative strand and the upstream a positive strand so we can flip the matrix
strand(CTCF_dFLAG_HMM) <- "*"
strand(CTCF_dFLAG_HMM[CTCF_dFLAG_HMM$updown == "upstream"]) <- "+"
strand(CTCF_dFLAG_HMM[CTCF_dFLAG_HMM$updown == "downstream"]) <- "-"



expNames <- ChIPdf[
    ChIPdf$cell %in% c("C17") &
        ChIPdf$prot %in% c("FLAG", "CTCF", "MAU2", "NIPBL", "SMC1") &
        ChIPdf$condition %in% c("TACL-ON", "TACL-OFF")
    # & ChIPdf$rep %in% paste0("rep", 1:10)
    , "name"
]

# > expNames
#  [1] "C17_TACL-ON_CTCF_rep1"       "C17_TACL-ON_FLAG_rep1"       "C17_TACL-ON_MAU2_rep1"       "C17_TACL-ON_FLAG_rep2"       "C17_TACL-ON_MAU2_rep2"
#  [6] "C17_TACL-ON_NIPBL_rep2"      "C17_TACL-ON_SMC1_rep2"       "C17_TACL-ON_SMC1_rep3"       "C17_TACL-ON_FLAG_rep1_rep2"  "C17_TACL-ON_MAU2_rep1_rep2"
# [11] "C17_TACL-ON_SMC1_rep2_rep3"  "C17_TACL-OFF_FLAG_rep2"      "C17_TACL-OFF_MAU2_rep2"      "C17_TACL-OFF_NIPBL_rep2"     "C17_TACL-OFF_SMC1_rep2"
# [16] "C17_TACL-OFF_SMC1_rep3"      "C17_TACL-OFF_SMC1_rep2_rep3"



ChIPdf_OI <- ChIPdf[ChIPdf$name %in% expNames, ]
any(file.exists(file.path(prefix, ChIPdf_OI$bigwig)))





toRnadoCov <- GetTornadoCov_V4(
    region = CTCF_dFLAG_HMM,
    ChIPdf = ChIPdf_OI,
    orderCov = c("C17_TACL-ON_FLAG_rep1_rep2", "C17_TACL-ON_MAU2_rep1_rep2"),
    orderby = "C17_TACL-ON_CTCF_rep1",
    zoom = 5000, binWidth = 10, Matrix = FALSE, makeList = FALSE, bincoverage = TRUE, prefix = prefix, bigwig = "pval"
)

saveRDS(toRnadoCov, file = "/storage/shared/TACL/PK/manuscript/revision/cov/Fig2_2E_CTCF_pval_240620.rds")





# first scale the data


toRnadoCov <- readRDS(file.path(covF, "Fig2_2E_CTCF_pval_240620.rds"))



tiles <- toRnadoCov$tiles
ChIPdf <- toRnadoCov$ChIPdf


for (i in seq_len(nrow(ChIPdf))) {
    name <- ChIPdf$name[i]
    meanBG <- ChIPdf$meanBG[i]
    cov_col <- paste0(name, "_cov")

    if (cov_col %in% colnames(mcols(tiles))) {
        # Scale the coverage values by dividing by meanBG
        message(paste("Scaling", cov_col, "by", meanBG))
        mcols(tiles)[[cov_col]] <- mcols(tiles)[[cov_col]] / meanBG
    }
}

toRnadoCov$tiles <- tiles

# Check the first few rows of the updated tiles
head(toRnadoCov$tiles)





# Add deltaCov
toRnadoCov$tiles$FLAG_cov <- toRnadoCov$tiles$"C17_TACL-ON_FLAG_rep1_rep2_cov" - toRnadoCov$tiles$"C17_TACL-OFF_FLAG_rep2_cov"
toRnadoCov$tiles$dMAU2_cov <- toRnadoCov$tiles$"C17_TACL-ON_MAU2_rep1_rep2_cov" - toRnadoCov$tiles$"C17_TACL-OFF_MAU2_rep2_cov"
toRnadoCov$tiles$dNIPBL_cov <- toRnadoCov$tiles$"C17_TACL-ON_NIPBL_rep2_cov" - toRnadoCov$tiles$"C17_TACL-OFF_NIPBL_rep2_cov"
toRnadoCov$tiles$dSMC1_cov <- toRnadoCov$tiles$"C17_TACL-ON_SMC1_rep2_rep3_cov" - toRnadoCov$tiles$"C17_TACL-OFF_SMC1_rep2_rep3_cov"




# add deltaCov to toRnadoCov$ChIPdf

# > colnames(toRnadoCov$ChIPdf)
#  [1] "name"            "manuscript"      "mapped_by"       "Lane"            "expname"         "RH_cell"         "cell"            "condition"
#  [9] "prot"            "rep"             "bigwig"          "fc_bw"           "meanBG"          "outerBG"         "TetO_enrichment"
# >



delta_ChIPdf <- data.frame(
    name = c("FLAG", "dMAU2", "dNIPBL", "dSMC1"),
    manuscript = rep(TRUE, 4),
    mapped_by = rep("PK", 4),
    Lane = rep(NA, 4),
    expname = c("dFLAG", "dMAU2", "dNIPBL", "dSMC1"),
    RH_cell = rep("C17_TMAU2", 4),
    cell = rep("C17", 4),
    condition = rep("ON-OFF", 4),
    prot = c("dFLAG", "dMAU2", "dNIPBL", "dSMC1"),
    rep = c("rep1_rep2", "rep1_rep2", "rep2", "rep2_rep3"),
    bigwig = rep(NA, 4),
    fc_bw = rep(NA, 4),
    meanBG = rep(1, 4),
    outerBG = rep(1, 4),
    TetO_enrichment = rep(NA, 4)
)

toRnadoCov$ChIPdf <- rbind(toRnadoCov$ChIPdf, delta_ChIPdf)








toRnadoCov$peaks$group <- NA

toRnadoCov$peaks[which(toRnadoCov$peaks$FIMO_convergent)]$group <- "HMM CTCF+;dFLAG+;convergent"

toRnadoCov$peaks[which(!toRnadoCov$peaks$FIMO_convergent)]$group <- "HMM CTCF+;dFLAG+;divergent"




GROUPSdf <- as.data.frame(table(toRnadoCov$peaks$group))
plotGroups <- GROUPSdf$Var1

# groups that will be plotted
GROUPSdf[GROUPSdf$Var1 %in% plotGroups, ]

# which groups will not be plotted
GROUPSdf[!GROUPSdf$Var1 %in% plotGroups, ]

# order bwOI per protein
expNames <- toRnadoCov$ChIPdf$name
bwOI_df <- data.frame(name = expNames, prot = toRnadoCov$ChIPdf[match(expNames, toRnadoCov$ChIPdf$name), "prot"], condition = toRnadoCov$ChIPdf[match(expNames, toRnadoCov$ChIPdf$name), "condition"])
prot_order <- c("CTCF", "FLAG", "dFLAG", "MAU2", "dMAU2", "NIPBL", "dNIPBL", "SMC1", "dSMC1")
bwOI_df$prot <- factor(bwOI_df$prot, levels = prot_order)
condition_order <- c("TACL-ON", "TACL-OFF", "ON-OFF")
bwOI_df$condition <- factor(bwOI_df$condition, levels = condition_order)
bwOI_df_ordered <- bwOI_df[order(bwOI_df$prot, bwOI_df$condition), ]
plotOrder <- bwOI_df_ordered$name
any(plotOrder %in% toRnadoCov$ChIPdf$name)

plotOrder <- c(
    bwOI_df_ordered$name
)

# order bwOI per protein
any(plotOrder %in% toRnadoCov$ChIPdf$name)




# matrix prenormalized before plotting

makePlotToRnado_V8(
    Plotregion_cov = toRnadoCov,
    plotOrder = plotOrder,
    orderFactor = "C17_TACL-ON_FLAG_rep1_rep2",
    plotGroups = plotGroups,
    norm = NULL,
    # norm = NULL,
    # normGroup = plotGroups,
    cutoffGroups = plotGroups,
    PlotMaxasCutoff = TRUE,
    SetCutOff = NULL,
    cutoffRatio = 1,
    cutoffLine = TRUE,
    subSampleGW = FALSE,
    # subSampleGroups = "plotGroups",
    # subsampleN=100,
    PeaksignalValue = 35,
    equalPlotMax = FALSE,
    setPlotMax = NULL,
    FlipRvStrand = TRUE,
    BigMatrix = TRUE,
    middleLine = TRUE,
    plotWidth = 500,
    protColors = protCols,
    PNG = FALSE,
    PDF = TRUE,
    outF = file.path(outF, "2E/"),
    plotName = paste0("2E_CTCF_dFLAG_HMM_signalVal35_", toRnadoCov$infoDF$bigwig, "_tornado_delta_noNorm"),
    ColWidth = 4
)

