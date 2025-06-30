#Fig 5

# Load packages -----------------------------------------------------------

library("rtracklayer")
library("GenomicRanges")

library("reshape2")
library("ggpubr")
library("dplyr")
library("ggplot2")
library("patchwork")

library(foreach)
library(doParallel)



# functions ----------------------------------------------------------------

# Add all functions to source file

prefix <- ""
source(file.path(prefix, "./TACL_functions_NG2025.r"))



# Features ----------------------------------------------------------------

# load TetO site locations
    # C17TetO <- read.table(file.path(prefix,'/TACL/PK/C17/final_C17locs_n27.bed'))
    # colnames(C17TetO) <- c('seqnames','start','end','score4C','judgement','TetO')
    # C17TetO <- GRanges(C17TetO)
    # saveRDS(C17TetO, file = '/storage/shared/TACL/PK/manuscript/revision/meta/final_C17locs_n27.rds')
C17TetO <- readRDS(file.path(prefix, "/TACL/PK/manuscript/revision/meta/final_C17locs_n27.rds"))

protCols <- read.table(file.path(prefix, "TACL/PK/manuscript/revision/meta/protCol_manuscript.tsv"), header = TRUE, comment.char = "")

dCols <- protCols[protCols$prot %in% c("FLAG", "MAU2", "SMC1", "RAD21"), ]
dCols$protein <- paste0("d", dCols$protein)
rbind(protCols, dCols)



# ChIP info
ChIPdf <- read.table(file.path(prefix, "TACL/PK/manuscript/revision/meta/ChIPdf_scalefactors_20250224.tsv"), header = TRUE, sep = "\t")

# Select only the ChIPs from the manuscript and remove the Input
ChIPdf <- ChIPdf[ChIPdf$manuscript, ]
ChIPdf <- ChIPdf[!ChIPdf$prot == "Input", ]


condition_levels <- c(
    "TACL-ON", "TACL-OFF", "mCherry", "WT","WT_IAA"
    "TACL-ON_IAA", "TACL-OFF_IAA", "mCherry_IAA",
    "TACL-ON_CRISPRi", "TACL-ON_CRISPRi-Ctrl", "TACL-ON_CRISPRi-WAPL",
    "TACL-ON_DRB", "mCherry_DRB"
)
rep_levels <- c("rep1", "rep2", "rep3", "rep4", "rep1_rep2", "rep1_rep3", "rep2_rep3", "rep1_rep2_rep3")

ChIPdf$condition <- factor(ChIPdf$condition, levels = condition_levels)
ChIPdf$rep <- factor(ChIPdf$rep, levels = rep_levels)
ChIPdf <- ChIPdf[order(ChIPdf$cell, ChIPdf$condition, ChIPdf$rep), ]



# add TetO enrichment
TetO_enrichment <- read.table(file.path(prefix, "TACL/PK/manuscript/revision/meta/TetO_enrichment_250224.tsv"), header = TRUE, sep = "\t")
ChIPdf <- merge(ChIPdf, TetO_enrichment[, c("name", "TetO_enrichment")], by = "name", all.x = TRUE)
ChIPdf <- ChIPdf[order(ChIPdf$cell, ChIPdf$condition, ChIPdf$rep), ]


# check whether all bigwigs exist
if (!all(file.exists(file.path(prefix, ChIPdf$bigwig)))) {
    message("Not all pval bigwigs exist")
    ChIPdf[!file.exists(file.path(prefix, ChIPdf$bigwig)), "name"]
}

if (!all(file.exists(file.path(prefix, ChIPdf$fc_bw)))) {
    message("Not all fc bigwigs exist")
    ChIPdf[!file.exists(file.path(prefix, ChIPdf$fc_bw)), "name"]
}



# HMM domains
domain_HMM <- read.table(file.path(prefix, "TACL/PK/manuscript/revision/meta/TACL_domains_HMM.bed"), header = FALSE, sep = "\t")
colnames(domain_HMM) <- c("chr", "start", "end", "name")
domain_HMM_GR <- GRanges(seqnames = domain_HMM$chr, ranges = IRanges(start = domain_HMM$start, end = domain_HMM$end))
# domain_HMM_GR$ID<-domain_HMM$name  #note mikhail added the wrong name for 2 and 25
domain_HMM_GR <- reduce(domain_HMM_GR) # remove redundant domains


# CTCF
peak_CTCF_C17_TMAU2 <- readRDS(file.path(prefix, "TACL/PK/CTCF_motif/C17_TACL_CTCF_hg38_peaks_withMotif_TetO.rds"))
peak_CTCF_C17_TMAU2$HMMdomain <- countOverlaps(peak_CTCF_C17_TMAU2, domain_HMM_GR)

peak_CTCF_C17_Tmcherry <- readRDS(file.path(prefix, "TACL/PK/CTCF_motif/C17_mCherry_CTCF_hg38_peaks_withMotif_TetO.rds"))


CTCF_reduced <- reduce(c(peak_CTCF_C17_TMAU2, peak_CTCF_C17_Tmcherry))

CTCF_reduced_sV35 <- reduce(c(peak_CTCF_C17_TMAU2[peak_CTCF_C17_TMAU2$signalValue > 35], peak_CTCF_C17_Tmcherry[peak_CTCF_C17_Tmcherry$signalValue > 35]))



peak_CTCF_C17_TMAU2_convergent <- peak_CTCF_C17_TMAU2[which(peak_CTCF_C17_TMAU2$FIMO_convergent)] # 14303
peak_CTCF_C17_TMAU2_divergent <- peak_CTCF_C17_TMAU2[which(peak_CTCF_C17_TMAU2$FIMO_convergent == FALSE)] # 14210



# C17 Flag peaks in TACL-ON
flag_opt <- import(file.path(prefix, "/TACL/PK/4DN_ChIP/outputs_all/chip.peaks/TACL_RH_eHAP1_C17_T-hMAU2v2_FLAG_mrgd-VER6376_rep1-VER7254_rep2_input_VER6376_idr.optimal_peak.regionPeak.bb")) # 11757
flag_opt_mrgdbyhighestPeak <- filterOverlapPeaks(peaks = flag_opt) # 10442
flag_opt_mrgdbyhighestPeak$ID <- 1:length(flag_opt_mrgdbyhighestPeak)
flag_opt_mrgdbyhighestPeak <- addTetOdistance(flag_opt_mrgdbyhighestPeak, C17TetO)
flag_opt_mrgdbyhighestPeak$HMMdomain <- countOverlaps(flag_opt_mrgdbyhighestPeak, domain_HMM_GR)

FLAG_filtered <- flag_opt_mrgdbyhighestPeak[flag_opt_mrgdbyhighestPeak$signalValue > 35] # 5725 peaks
table(FLAG_filtered$distance < 3e6) # F 4795; T 930

# dFLAG peaks
dFLAG <- readRDS(file.path(prefix, "TACL/PK/manuscript/revision/meta/dFLAG_peaks.rds"))

# NKI peaks
NKI_H3K27ac_peaks <- import(file.path(prefix, "TACL/NKI/peaks/all_H3K27ac_peaks.canonical.replicated.no_blacklist.merged.bed"))
NKI_H3K4me3_peaks <- import(file.path(prefix, "TACL/NKI/peaks/all_H3K4me3_peaks.canonical.replicated.no_blacklist.merged.bed"))

# H3K4me1 (>=1)
H3K4me1_opt_All <- import(paste0(prefix, "TACL/PK/Public_ChIP/Encode_H3K4me1/HAP1_H3K4me1_ENCODE_hg38_pseudorep_ENCFF963OJW.bb"))
H3K4me1_opt_All <- H3K4me1_opt_All[H3K4me1_opt_All$signalValue >= 1]
H3K4me1_opt_reduced <- GenomicRanges::reduce(H3K4me1_opt_All) # 116897



# STAG2 lost domains
STAG2_lost_domains <- read.table("/storage/shared/TACL/PK/manuscript/revision/meta/STAG2_lost_domains.txt", header = TRUE, sep = "\t")
STAG2_lost_domains_GR <- GRanges(seqnames = STAG2_lost_domains$chrom, ranges = IRanges(start = STAG2_lost_domains$tacl_domain_stag2_iaa_start, end = STAG2_lost_domains$tacl_domain_stag2_iaa_end))
STAG2_lost_domains_GR$teto_id <- STAG2_lost_domains$teto_id

STAG2_lost_domains_WT_GR <- GRanges(seqnames = STAG2_lost_domains$chrom, ranges = IRanges(start = STAG2_lost_domains$tacl_domain_stag2_wt_start, end = STAG2_lost_domains$tacl_domain_stag2_wt_end))
STAG2_lost_domains_WT_GR$teto_id <- STAG2_lost_domains$teto_id




# degron domains
degron_domains <- read.table("/storage/shared/TACL/PK/manuscript/revision/meta/degrons_domains_HMM.txt", header = TRUE, sep = "\t")
#chrom	start	end	teto_id	tacl_domain_ctcf_start	tacl_domain_ctcf_end	tacl_domain_wapl_start	tacl_domain_wapl_end	tacl_domain_stag2_start	tacl_domain_stag2_end	tacl_domain_pds5a_start	tacl_domain_pds5a_end

CTCF_degron_domain_GR <- GRanges(seqnames = degron_domains$chrom, ranges = IRanges(start = degron_domains$tacl_domain_ctcf_start, end = degron_domains$tacl_domain_ctcf_end))
CTCF_degron_domain_GR$ID <- degron_domains$teto_id

degron_domains_WAPL <- degron_domains[!is.na(degron_domains$tacl_domain_wapl_start), ]
WAPL_degron_domain_GR <- GRanges(seqnames = degron_domains_WAPL$chrom, ranges = IRanges(start = degron_domains_WAPL$tacl_domain_wapl_start, end = degron_domains_WAPL$tacl_domain_wapl_end))
WAPL_degron_domain_GR$ID <- degron_domains_WAPL$teto_id

STAG2_degron_domain_GR <- GRanges(seqnames = degron_domains$chrom, ranges = IRanges(start = degron_domains$tacl_domain_stag2_start, end = degron_domains$tacl_domain_stag2_end))
STAG2_degron_domain_GR$ID <- degron_domains$teto_id

PDS5A_degron_domain_GR <- GRanges(seqnames = degron_domains$chrom, ranges = IRanges(start = degron_domains$tacl_domain_pds5a_start, end = degron_domains$tacl_domain_pds5a_end))
PDS5A_degron_domain_GR$ID <- degron_domains$teto_id





#4C

normalized_4C_rep2 <- readRDS("/storage/shared/TACL/PK/4C/TACL_rep2_normalized_4C.rds")
normalized_4C <- readRDS("/storage/shared/TACL/PK/4C/TACL_C17_degrons_including_mch_VER10253_normalized_4C.rds")
normalized_4C_RAD21 <- readRDS("/storage/shared/TACL/PK/4C/TACL_C17_RAD21AID_normalized_4C_240709.rds")





# #CTCF peaks ---------------------------------------------------------------

# Give the downstream TetO motifs a negative strand and the upstream a positive strand so we can flip the matrix
strand(peak_CTCF_C17_TMAU2) <- "*"
strand(peak_CTCF_C17_TMAU2[peak_CTCF_C17_TMAU2$updown == "upstream"]) <- "+"
strand(peak_CTCF_C17_TMAU2[peak_CTCF_C17_TMAU2$updown == "downstream"]) <- "-"



#intersect CTCF peaks with STAG2 lost domains

peak_CTCF_C17_TMAU2$STAG2_lost_domain <- countOverlaps(peak_CTCF_C17_TMAU2, STAG2_lost_domains_GR)


#Add extended degron domains
peak_CTCF_C17_TMAU2$degron_domain_CTCF <- countOverlaps(peak_CTCF_C17_TMAU2, CTCF_degron_domain_GR)
peak_CTCF_C17_TMAU2$degron_domain_WAPL <- countOverlaps(peak_CTCF_C17_TMAU2, WAPL_degron_domain_GR)
peak_CTCF_C17_TMAU2$degron_domain_STAG2 <- countOverlaps(peak_CTCF_C17_TMAU2, STAG2_degron_domain_GR)
peak_CTCF_C17_TMAU2$degron_domain_PDS5A <- countOverlaps(peak_CTCF_C17_TMAU2, PDS5A_degron_domain_GR)



#ChIP coverage
peak_CTCF_C17_TMAU2<-getCov(bw=file.path(prefix,ChIPdf[ChIPdf$name == "C17_TACL-ON_CTCF_rep1", "bigwig"]), region=peak_CTCF_C17_TMAU2, name="C17_TACL-ON_CTCF_rep1_cov")
peak_CTCF_C17_TMAU2<-getCov(bw=file.path(prefix,ChIPdf[ChIPdf$name == "C17_TACL-ON_FLAG_rep1_rep2", "bigwig"]), region=peak_CTCF_C17_TMAU2, name="C17_TACL-ON_FLAG_rep1_rep2_cov")
peak_CTCF_C17_TMAU2<-getCov(bw=file.path(prefix,ChIPdf[ChIPdf$name == "C17_TACL-ON_MAU2_rep1_rep2", "bigwig"]), region=peak_CTCF_C17_TMAU2, name="C17_TACL-ON_MAU2_rep1_rep2_cov")
peak_CTCF_C17_TMAU2<-getCov(bw=file.path(prefix,ChIPdf[ChIPdf$name == "C17_TACL-ON_NIPBL_rep2", "bigwig"]), region=peak_CTCF_C17_TMAU2, name="C17_TACL-ON_NIPBL_rep2_cov")
peak_CTCF_C17_TMAU2<-getCov(bw=file.path(prefix,ChIPdf[ChIPdf$name == "C17_TACL-ON_RAD21_rep1_rep3", "bigwig"]), region=peak_CTCF_C17_TMAU2, name="C17_TACL-ON_RAD21_rep1_rep3_cov")
peak_CTCF_C17_TMAU2<-getCov(bw=file.path(prefix,ChIPdf[ChIPdf$name == "C17_TACL-ON_SMC1_rep2_rep3", "bigwig"]), region=peak_CTCF_C17_TMAU2, name="C17_TACL-ON_SMC1_rep2_rep3_cov")
peak_CTCF_C17_TMAU2<-getCov(bw=file.path(prefix,ChIPdf[ChIPdf$name == "C17_TACL-ON_STAG1_rep1", "bigwig"]), region=peak_CTCF_C17_TMAU2, name="C17_TACL-ON_STAG1_rep1_cov")
peak_CTCF_C17_TMAU2<-getCov(bw=file.path(prefix,ChIPdf[ChIPdf$name == "C17_TACL-ON_STAG2_rep1", "bigwig"]), region=peak_CTCF_C17_TMAU2, name="C17_TACL-ON_STAG2_rep1_cov")
peak_CTCF_C17_TMAU2<-getCov(bw=file.path(prefix,ChIPdf[ChIPdf$name == "C17_WT_STAG1_rep1", "bigwig"]), region=peak_CTCF_C17_TMAU2, name="C17_WT_STAG1_rep1_cov")




#group based on TACl domain and STAG2_lost_domain



mcols(peak_CTCF_C17_TMAU2)$groupHMM <- case_when(
  mcols(peak_CTCF_C17_TMAU2)$HMMdomain == 0 & mcols(peak_CTCF_C17_TMAU2)$distance > 3e6 ~ "GW",
  mcols(peak_CTCF_C17_TMAU2)$HMMdomain > 0 ~ "HMM",
  TRUE ~ NA_character_  # Set NA for all other cases
)

# Step 1: Add "group" column to the GRanges metadata
mcols(peak_CTCF_C17_TMAU2)$group <- case_when(
  mcols(peak_CTCF_C17_TMAU2)$HMMdomain == 0 & mcols(peak_CTCF_C17_TMAU2)$distance > 3e6 ~ "GW",
  mcols(peak_CTCF_C17_TMAU2)$HMMdomain > 0 ~ "HMM",
  TRUE ~ NA_character_  # Set NA for all other cases
)

# Step 2: Extract unique "HMM" groups from the peaks
groups <- unique(mcols(peak_CTCF_C17_TMAU2)$group)
groups <- groups[!is.na(groups)]
groups <- groups[grep("HMM", groups)]  # Select only HMM-related groups

# Step 3: Update the "group" column based on "STAG2_lost_domain"
for (group in groups) {
  message(group)  # Print the group being processed

  # Find indices where group is "HMM" and STAG2_lost_domain > 0
  IDX <- which(mcols(peak_CTCF_C17_TMAU2)$group %in% group & mcols(peak_CTCF_C17_TMAU2)$STAG2_lost_domain > 0)
  if (length(IDX) > 0) {
    mcols(peak_CTCF_C17_TMAU2)$group[IDX] <- paste0(group, ";STAG2 lost domain+")
  }

  # Find indices where group is "HMM" and STAG2_lost_domain == 0
  IDX <- which(mcols(peak_CTCF_C17_TMAU2)$group %in% group & mcols(peak_CTCF_C17_TMAU2)$STAG2_lost_domain == 0)
  if (length(IDX) > 0) {
    mcols(peak_CTCF_C17_TMAU2)$group[IDX] <- paste0(group, ";STAG2 lost domain-")
  }
}

# Step 4: Check the updated metadata
mcols(peak_CTCF_C17_TMAU2)






#Within these domains, which we distinguished three categories of CTCF sites with distinct ChIP-seq signals: the strong, intermediate and weak CTCF binding sites (Fig. xx).



# Step 1: Compute quantiles (33% and 66%) for all CTCF peaks
quantiles_all <- quantile(peak_CTCF_C17_TMAU2$'C17_TACL-ON_CTCF_rep1_cov', probs = c(0.33, 0.66), na.rm = TRUE)

# Step 2: Compute quantiles (33% and 66%) for HMM-specific CTCF peaks
HMM_IDX <- which(mcols(peak_CTCF_C17_TMAU2)$HMMdomain > 0)

if (length(HMM_IDX) > 0) {
  quantiles_HMM <- quantile(mcols(peak_CTCF_C17_TMAU2)$'C17_TACL-ON_CTCF_rep1_cov'[HMM_IDX], 
                             probs = c(0.33, 0.66), 
                             na.rm = TRUE)
} else {
  message("No peaks found where HMMdomain > 0")
  quantiles_HMM <- c(NA, NA)  # Assign NA if no HMM peaks exist
}

# Step 3: Add new columns for classification
mcols(peak_CTCF_C17_TMAU2)$CTCF_strength_all <- case_when(
  mcols(peak_CTCF_C17_TMAU2)$'C17_TACL-ON_CTCF_rep1_cov' <= quantiles_all[1] ~ "Low",
  mcols(peak_CTCF_C17_TMAU2)$'C17_TACL-ON_CTCF_rep1_cov' > quantiles_all[1] & 
    mcols(peak_CTCF_C17_TMAU2)$'C17_TACL-ON_CTCF_rep1_cov' <= quantiles_all[2] ~ "Medium",
  mcols(peak_CTCF_C17_TMAU2)$'C17_TACL-ON_CTCF_rep1_cov' > quantiles_all[2] ~ "High"
)

# Initialize column with NA
mcols(peak_CTCF_C17_TMAU2)$CTCF_strength_HMM <- NA_character_

if (length(HMM_IDX) > 0) {
  mcols(peak_CTCF_C17_TMAU2)$CTCF_strength_HMM <- case_when(
    mcols(peak_CTCF_C17_TMAU2)$'C17_TACL-ON_CTCF_rep1_cov' <= quantiles_HMM[1] ~ "Low",
    mcols(peak_CTCF_C17_TMAU2)$'C17_TACL-ON_CTCF_rep1_cov' > quantiles_HMM[1] & 
      mcols(peak_CTCF_C17_TMAU2)$'C17_TACL-ON_CTCF_rep1_cov' <= quantiles_HMM[2] ~ "Medium",
    mcols(peak_CTCF_C17_TMAU2)$'C17_TACL-ON_CTCF_rep1_cov' > quantiles_HMM[2] ~ "High"
  )
}


table(peak_CTCF_C17_TMAU2$CTCF_strength_all, peak_CTCF_C17_TMAU2$CTCF_strength_HMM)
  #        High   Low Medium
  # High   12369     0      0
  # Low        0 12006      0
  # Medium   490    59  11456



saveRDS(peak_CTCF_C17_TMAU2, file = "/storage/shared/TACL/PK/manuscript/revision/fig2_revision/CTCFpeaks_degrons/CTCF_peaks_CTCF_250327.rds")








#ChIP coverage in/out collapsed domains and GW

#plot coverage distribution
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

peak_CTCF_C17_TMAU2 <- readRDS("/storage/shared/TACL/PK/manuscript/revision/fig2_revision/CTCFpeaks_degrons/CTCF_peaks_CTCF_250327.rds")


table(peak_CTCF_C17_TMAU2$signalValue>35, peak_CTCF_C17_TMAU2$HMMdomain>0)
#       FALSE  TRUE
#  FALSE  8260   179
#  TRUE  27329   612   > 35



#boxplots for CTCF coverage in collapsed domains and GW -----------------------------------------------------

# Step 1: Convert GRanges metadata to a data frame
plot_data <- mcols(peak_CTCF_C17_TMAU2) %>%
  as.data.frame() %>%
  select(
    group, CTCF_strength_HMM,
    C17_TACL.ON_CTCF_rep1_cov, C17_TACL.ON_FLAG_rep1_rep2_cov, 
    C17_TACL.ON_MAU2_rep1_rep2_cov, C17_TACL.ON_NIPBL_rep2_cov,
    C17_TACL.ON_RAD21_rep1_rep3_cov, C17_TACL.ON_SMC1_rep2_rep3_cov,
    C17_TACL.ON_STAG1_rep1_cov, C17_TACL.ON_STAG2_rep1_cov, C17_WT_STAG1_rep1_cov
  ) %>%
  pivot_longer(cols = matches("^C17_TACL\\.ON|^C17_WT"),  # Corrected regex to match both prefixes
               names_to = "Coverage_Type", values_to = "Coverage_Value") %>%
  filter(!is.na(group))  # Exclude NA groups

# Step 2: Define color scheme for CTCF strength
ctcf_colors <- c("Low" = "blue", "Medium" = "orange", "High" = "red")

# Step 3: Generate and save boxplots for each coverage type
coverage_types <- unique(plot_data$Coverage_Type)

for (coverage in coverage_types) {
  
message(coverage)  # Print the coverage type being processed

  # Filter data for the current coverage type
  coverage_data <- plot_data %>% filter(Coverage_Type == coverage)

  coverage_data$CTCF_strength_HMM <- factor(
  coverage_data$CTCF_strength_HMM,
  levels = c("Low", "Medium", "High")
)

  
  #scale using meanBG
  protein_scalefactor <- sub("TACL.ON","TACL-ON",coverage)
  protein_scalefactor <- sub("_cov","",protein_scalefactor)
  scalefactor <- ChIPdf[ChIPdf$name == protein_scalefactor, "meanBG"]
  coverage_data$Coverage_meanBG <- coverage_data$Coverage_Value/scalefactor
  
# Calculate dynamic y-axis limits (ignore extreme outliers)
y_min <- min(coverage_data$Coverage_meanBG, na.rm = TRUE)
#y_max <- quantile(coverage_data$Coverage_meanBG, probs = 0.99, na.rm = TRUE)  # Use the 99th percentile

# Calculate Q3 and IQR
# Compute whisker max (Q3 + 1.5*IQR) per group and strength
library(dplyr)

group_stats <- coverage_data %>%
  group_by(group, CTCF_strength_HMM) %>%
  summarise(
    Q3 = quantile(Coverage_meanBG, 0.75, na.rm = TRUE),
    Q1 = quantile(Coverage_meanBG, 0.25, na.rm = TRUE),
    IQR = Q3 - Q1,
    whisker_max = Q3 + 1.5 * IQR,
    .groups = "drop"
  )

# Global y_max = largest whisker max across groups, with 10% buffer
y_max <- ceiling(max(group_stats$whisker_max, na.rm = TRUE))



boxplot_ctcf <- ggplot(coverage_data, aes(x = group, y = Coverage_meanBG, fill = CTCF_strength_HMM)) +
  geom_boxplot(outlier.shape = NA, position = "dodge", color = "black") +  # Remove outliers
  coord_cartesian(ylim = c(y_min, y_max)) +  # Auto-scale y-axis using the 99th percentile
  theme_minimal() +
  scale_fill_manual(values = ctcf_colors) +
  #geom_jitter(color="black", size=0.4, alpha=0.9) +
  labs(
    title = paste("CTCF Peak Coverage -", coverage),
    x = "Group",
    y = "Normalized Coverage Value (MeanBG)",
    fill = "CTCF Strength (HMM)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the boxplot
ggsave(paste0("/storage/shared/TACL/PK/manuscript/revision/fig2_revision/CTCFpeaks_degrons/ChIPcov/250327_collapsed_CTCF_Peak_Coverage_Boxplot_", coverage, ".pdf"), boxplot_ctcf, width = 8, height = 6)


}






#ASP plots ---------------------------------------


peak_CTCF_C17_TMAU2 <- readRDS("/storage/shared/TACL/PK/manuscript/revision/fig2_revision/CTCFpeaks_degrons/CTCF_peaks_CTCF_250327.rds")


expNames <- c(
    "C17_TACL-ON_FLAG_rep1_rep2", "C17_TACL-OFF_FLAG_rep2", "C17_mCherry_FLAG_rep1",
    "C17_TACL-ON_MAU2_rep1_rep2", "C17_TACL-OFF_MAU2_rep2", "C17_mCherry_MAU2_rep1_rep2",
    "C17_TACL-ON_NIPBL_rep2", "C17_TACL-OFF_NIPBL_rep2", "C17_mCherry_NIPBL_rep1_rep2",
    "C17_TACL-ON_RAD21_rep1_rep3", "C17_TACL-OFF_RAD21_rep3", "C17_mCherry_RAD21_rep1_rep3",
    "C17_TACL-ON_SMC1_rep2_rep3", "C17_TACL-OFF_SMC1_rep2_rep3", "C17_mCherry_SMC1_rep2_rep3",
    "C17_TACL-ON_STAG1_rep1", "C17_TACL-OFF_STAG1_rep1", "C17_mCherry_STAG1_rep1",
    "C17_TACL-ON_STAG2_rep1", "C17_TACL-OFF_STAG2_rep1", "C17_mCherry_STAG2_rep1",
    "C17_TACL-ON_CTCF_rep1", "C17_mCherry_CTCF_rep1",
    "C17_TACL-ON_PDS5A_rep1", "C17_mCherry_PDS5A_rep1",
    "C17_TACL-ON_WAPL_rep1", "C17_mCherry_WAPL_rep1",
    "C17_STAG2AID_C2_TACL-ON_IAA_NIPBL_rep1",
    "C17_STAG2AID_C2_TACL-ON_IAA_SMC1_rep1",
    "C17_STAG2AID_C2_TACL-ON_IAA_STAG1_rep1",
    "C17_STAG2AID_C2_TACL-ON_IAA_FLAG_rep1_rep2"
)

#TACL-ON-SA2-AID
#ChIPdf[ChIPdf$cell=="C17_STAG2AID_C2" & ChIPdf$condition=="TACL-ON_IAA",]

ChIPdf_OI <- ChIPdf[ChIPdf$name %in% expNames, ]



toRnadoCov <- GetTornadoCov_V4(
  region = peak_CTCF_C17_TMAU2,
  ChIPdf = ChIPdf_OI,
  orderCov = "C17_TACL-ON_CTCF_rep1",
  orderby = "C17_TACL-ON_CTCF_rep1",
  zoom = 5000, binWidth = 10, Matrix = FALSE, makeList = FALSE, bincoverage = TRUE, prefix = prefix, bigwig = "pval"
)

saveRDS(toRnadoCov, file = "/storage/shared/TACL/PK/manuscript/revision/cov/Fig5_CTCFpeaks_pval_250331.rds")


toRnadoCov <- readRDS("/storage/shared/TACL/PK/manuscript/revision/cov/Fig5_CTCFpeaks_pval_250331.rds")

toRnadoCov$peaks$groupV2 <- toRnadoCov$peaks$group

toRnadoCov$peaks$convergent <- ifelse(
  is.na(toRnadoCov$peaks$FIMO_convergent), 
  NA, 
  ifelse(toRnadoCov$peaks$FIMO_convergent, "convergent", "divergent")
)


#unique(toRnadoCov$peaks$group)
#[1] "GW"                     NA                       "HMM;STAG2 lost domain+" "HMM;STAG2 lost domain-"


toRnadoCov$peaks$group <- paste0(toRnadoCov$peaks$groupV2, ";", toRnadoCov$peaks$CTCF_strength_all,";", toRnadoCov$peaks$convergent)

unique(toRnadoCov$peaks$group)
# [1] "GW;High"                       "NA;High"                       "HMM;STAG2 lost domain+;High"   "HMM;STAG2 lost domain-;High"
# [5] "GW;Medium"                     "NA;Medium"                     "HMM;STAG2 lost domain-;Medium" "HMM;STAG2 lost domain+;Medium"
# [9] "GW;Low"                        "NA;Low"                        "HMM;STAG2 lost domain-;Low"    "HMM;STAG2 lost domain+;Low"

plotGroups <- unique(toRnadoCov$peaks$group)[!grepl("^NA;|;NA", unique(toRnadoCov$peaks$group))]


table(toRnadoCov$peaks[toRnadoCov$peaks$group %in% plotGroups, ]$group, toRnadoCov$peaks[toRnadoCov$peaks$group %in% plotGroups, ]$CTCF_strength_HMM)
$group)



plotZoom <- 5000
plotWidth <- 1000
binWidth <- 10
# V151-V350
avgdf <- avgPlot_V2(toRnadoCov,
    plotOrder = NULL,
    #plotOrder = plotOrder[grepl(pattern = paste(ChIP_factors, collapse = "|"), plotOrder)],
    plotGroups = plotGroups,
    orderFactor = NULL,
    norm = "meanBG",
    #normGroup = plotGroups,
    PeaksignalValue = 0,
    plotWidth = plotWidth,
    FlipRvStrand = TRUE,
    middleLine = TRUE
)

library(dplyr)

avgdf <- avgdf %>%
  mutate(
    grouptype = case_when(
      grepl("^GW;", group) ~ "GW",
      grepl("^HMM;STAG2 lost domain\\+;", group) ~ "HMM;STAG2 lost domain+",
      grepl("^HMM;STAG2 lost domain-;", group) ~ "HMM;STAG2 lost domain-",
      TRUE ~ NA_character_
    ),
    level = case_when(
      grepl(";High;", group) ~ "High",
      grepl(";Medium;", group) ~ "Medium",
      grepl(";Low;", group) ~ "Low",
      TRUE ~ NA_character_
    ),
    orientation = case_when(
      grepl("convergent$", group) ~ "convergent",
      grepl("divergent$", group) ~ "divergent",
      TRUE ~ NA_character_
    )
  )

# Optional: convert to factors for plotting
avgdf$level <- factor(avgdf$level, levels = c("Low", "Medium", "High"))
avgdf$orientation <- factor(avgdf$orientation, levels = c("convergent", "divergent"))


unique(avgdf$group)
unique(avgdf$grouptype)
unique(avgdf$level)
unique(avgdf$orientation)


middleBins <- c(plotZoom / binWidth / 2, (plotZoom / binWidth / 2) + 1)


#library(tidyverse)
library(dplyr)
library(tidyr)      # for `separate()`, `pivot_longer()`
library(ggplot2)    # for plotting
library(patchwork)  # for combining plots

outF<-"/storage/shared/TACL/PK/manuscript/revision/fig2_revision/CTCFpeaks_degrons/ASP/"

exp_names <- unique(avgdf$expName)
grouptypes <- c("GW", "HMM;STAG2 lost domain+", "HMM;STAG2 lost domain-")
orientations <- c("convergent", "divergent")




#First save as list

plots_by_expname <- list()

for (exp in exp_names) {

  plot_list <- list()

  # Calculate global y-axis range for this expName
  global_y_range <- avgdf %>%
    filter(expName == exp) %>%
    pivot_longer(cols = starts_with("V"), names_to = "Position", values_to = "Value") %>%
    summarize(ymin = min(Value, na.rm = TRUE), ymax = max(Value, na.rm = TRUE))

  for (ori in orientations) {
    for (grp in grouptypes) {

      df_plot <- avgdf %>% filter(expName == exp, grouptype == grp, orientation == ori)

      if (nrow(df_plot) == 0) next

      df_long <- df_plot %>%
        pivot_longer(cols = starts_with("V"), names_to = "Position", values_to = "Value") %>%
        mutate(
          Position = as.numeric(gsub("V", "", Position)),
          New_Position = (Position - middleBins[1]) * binWidth
        )

      p <- ggplot(df_long, aes(x = New_Position, y = Value, color = level)) +
        geom_line(linewidth = 1.2) +
        labs(
          title = paste(grp, ori, sep = " | "),
          x = "Position (bp)",
          y = "Normalized Signal"
        ) +
        scale_color_manual(values = c("Low" = "#a6cee3", "Medium" = "#1f78b4", "High" = "#08306b")) +
        coord_cartesian(ylim = c(global_y_range$ymin, global_y_range$ymax)) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
          legend.position = "top",
          legend.title = element_blank()
        )

      plot_key <- paste(ori, grp, sep = "_")
      plot_list[[plot_key]] <- p
    }
  }

  
  
  # Save all subplots for this expName
  plots_by_expname[[exp]] <- plot_list

  # Combine selected subplots into one grid and save as PDF
  ordered_keys <- c(
    paste0("convergent_", grouptypes),
    paste0("divergent_", grouptypes)
  )

  combined_plot <- wrap_plots(plot_list[ordered_keys], ncol = 3, byrow = TRUE) +
    plot_annotation(title = exp) &
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

  ggsave(
    file.path(outF, paste0(gsub("[^A-Za-z0-9]", "_", exp), "_plots.pdf")),
    combined_plot,
    width = 16, height = 10
  )
}


# View a specific subplot
plots_by_expname[["C17_TACL-ON_FLAG_rep1_rep2"]][["convergent_HMM;STAG2 lost domain+"]]

# Print it
print(plots_by_expname[[1]][[1]])






##########
# Function to create overview grid of ASP plots

library(ggplot2)
library(patchwork)
library(stringr)
library(grid)




plot_overview_grid <- function(overview_exps, plots_by_expname, outF, outfile_name, title = "ASP Signal by Orientation and Domain Loss", y_axis_ranges = NULL, remove_titles = TRUE) {
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(grid)

  # Extract protein names from expName dynamically using improved regex
  protein_labels <- stringr::str_extract(basename(overview_exps), "(?<=_)(FLAG|NIPBL|SMC1|STAG1|STAG2)(?=_?rep|$)")

  # Define row labels (2-line format)
  row_keys <- c(
    "convergent_HMM;STAG2 lost domain+",
    "convergent_HMM;STAG2 lost domain-",
    "divergent_HMM;STAG2 lost domain+",
    "divergent_HMM;STAG2 lost domain-"
  )
  row_labels <- c(
    "convergent\ninside (+)",
    "convergent\noutside (-)",
    "divergent\ninside (+)",
    "divergent\noutside (-)"
  )

  # Build grid of plots from provided lookup
  overview_plots <- list()
  for (exp in overview_exps) {
    for (row_key in row_keys) {
      plot_key <- paste(exp, row_key, sep = "_")
      if (!is.null(plots_by_expname[[exp]][[row_key]])) {
        p <- plots_by_expname[[exp]][[row_key]]

        # Handle subplot titles
        if (remove_titles) {
          p <- p + labs(title = NULL) + theme(plot.title = element_blank())
        } else {
          p <- p + labs(title = row_key) + theme(plot.title = element_text(size = 8, hjust = 0.5))
        }

        # Apply y-axis range if specified
        if (!is.null(y_axis_ranges) && !is.null(y_axis_ranges[[exp]])) {
          range_vals <- y_axis_ranges[[exp]]
          p <- p + coord_cartesian(ylim = range_vals)
        }

        overview_plots[[plot_key]] <- p
      } else {
        overview_plots[[plot_key]] <- ggplot() + theme_void()
      }
    }
  }

  # Build labeled grid
  plot_grid <- list()
  for (i in seq_along(row_keys)) {
    row <- list()

    # Add row label (rotated on left)
    row_label <- wrap_elements(grid::textGrob(row_labels[i], rot = 90, gp = gpar(fontsize = 11, fontface = "bold")))
    row[[1]] <- row_label

    for (j in seq_along(overview_exps)) {
      key <- paste(overview_exps[j], row_keys[i], sep = "_")
      row[[j + 1]] <- overview_plots[[key]]
    }

    plot_grid[[i]] <- wrap_plots(row, ncol = length(row), widths = c(1.5, rep(5, length(overview_exps))))
  }

  # Create column headers (protein names)
  column_titles <- lapply(protein_labels, function(lbl) {
    ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = lbl, size = 6, fontface = "bold", hjust = 0.5, vjust = 0.5)
  })
  blank_plot <- ggplot() + theme_void()
  column_titles <- c(list(blank_plot), column_titles)
  column_header <- wrap_plots(column_titles, ncol = length(column_titles), widths = c(1.5, rep(5, length(overview_exps))))

  # Combine header, grid, and shared legend
  overview_plot <- column_header / wrap_plots(plot_grid, ncol = 1, guides = "collect") +
    plot_layout(heights = c(1, 20)) +
    plot_annotation(title = title) &
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )

  # Adjust output width based on number of columns (standard width per column)
  col_width <- 4
  width_total <- 1.5 + length(overview_exps) * col_width

  # Save
  ggsave(
    file.path(outF, outfile_name),
    overview_plot,
    width = width_total,
    height = 14
  )

  return(overview_plot)
}


> names(plots_by_expname)
 [1] "C17_mCherry_CTCF_rep1"                      "C17_mCherry_FLAG_rep1"                      "C17_mCherry_PDS5A_rep1"
 [4] "C17_mCherry_STAG1_rep1"                     "C17_mCherry_STAG2_rep1"                     "C17_mCherry_WAPL_rep1"
 [7] "C17_mCherry_MAU2_rep1_rep2"                 "C17_mCherry_NIPBL_rep1_rep2"                "C17_mCherry_RAD21_rep1_rep3"
[10] "C17_mCherry_SMC1_rep2_rep3"                 "C17_TACL-OFF_STAG1_rep1"                    "C17_TACL-OFF_STAG2_rep1"
[13] "C17_TACL-OFF_FLAG_rep2"                     "C17_TACL-OFF_MAU2_rep2"                     "C17_TACL-OFF_NIPBL_rep2"
[16] "C17_TACL-OFF_RAD21_rep3"                    "C17_TACL-OFF_SMC1_rep2_rep3"                "C17_TACL-ON_CTCF_rep1"
[19] "C17_TACL-ON_PDS5A_rep1"                     "C17_TACL-ON_STAG1_rep1"                     "C17_TACL-ON_STAG2_rep1"
[22] "C17_TACL-ON_WAPL_rep1"                      "C17_TACL-ON_NIPBL_rep2"                     "C17_TACL-ON_FLAG_rep1_rep2"
[25] "C17_TACL-ON_MAU2_rep1_rep2"                 "C17_TACL-ON_RAD21_rep1_rep3"                "C17_TACL-ON_SMC1_rep2_rep3"
[28] "C17_STAG2AID_C2_TACL-ON_IAA_NIPBL_rep1"     "C17_STAG2AID_C2_TACL-ON_IAA_SMC1_rep1"      "C17_STAG2AID_C2_TACL-ON_IAA_STAG1_rep1"
[31] "C17_STAG2AID_C2_TACL-ON_IAA_FLAG_rep1_rep2"



#TACL-ON

overview_exps <- c(
  "C17_TACL-ON_FLAG_rep1_rep2",
  "C17_TACL-ON_NIPBL_rep2",
  "C17_TACL-ON_SMC1_rep2_rep3",
  "C17_TACL-ON_STAG1_rep1",
  "C17_TACL-ON_STAG2_rep1"
)

y_axis_ranges <- list(
  "C17_TACL-ON_FLAG_rep1_rep2" = c(0, 7),
  "C17_TACL-ON_NIPBL_rep2" = c(0, 20),
  "C17_TACL-ON_SMC1_rep2_rep3" = c(0, 3),
  "C17_TACL-ON_STAG1_rep1" = c(0, 3),
  "C17_TACL-ON_STAG2_rep1" = c(0, 4)
)

TACLON_ASP <- plot_overview_grid(
  overview_exps = overview_exps,
  plots_by_expname = plots_by_expname,
  outF = outF,
  outfile_name = "250402_TACL-ON_ASP_overview_plot.pdf",
  title = "TACL-ON",
  y_axis_ranges = y_axis_ranges,
  remove_titles = TRUE
)

plot_overview_grid(
  overview_exps = overview_exps,
  plots_by_expname = plots_by_expname,
  outF = outF,
  outfile_name = "250402_TACL-ON_ASP_overview_plot_withTitle.pdf",
  title = "TACL-ON",
  y_axis_ranges = y_axis_ranges,
  remove_titles = FALSE
)


#TACL-OFF

overview_exps <- c(
  "C17_TACL-OFF_FLAG_rep2",
  "C17_TACL-OFF_NIPBL_rep2",
  "C17_TACL-OFF_SMC1_rep2_rep3",
  "C17_TACL-OFF_STAG1_rep1",
  "C17_TACL-OFF_STAG2_rep1"
)

y_axis_ranges <- list(
  "C17_TACL-OFF_FLAG_rep2" = c(0, 7),
  "C17_TACL-OFF_NIPBL_rep2" = c(0, 20),
  "C17_TACL-OFF_SMC1_rep2_rep3" = c(0, 3),
  "C17_TACL-OFF_STAG1_rep1" = c(0, 3),
  "C17_TACL-OFF_STAG2_rep1" = c(0, 4)
)

TACLOFF_ASP <- plot_overview_grid(
  overview_exps = overview_exps,
  plots_by_expname = plots_by_expname,
  outF = outF,
  outfile_name = "250409_TACL-OFF_ASP_overview_plot.pdf",
  title = "TACL-OFF",
  y_axis_ranges = y_axis_ranges,
  remove_titles = TRUE
)

plot_overview_grid(
  overview_exps = overview_exps,
  plots_by_expname = plots_by_expname,
  outF = outF,
  outfile_name = "250409_TACL-OFF_ASP_overview_plot_withTitle.pdf",
  title = "TACL-OFF",
  y_axis_ranges = y_axis_ranges,
  remove_titles = FALSE
)

#mCherry

overview_exps <- c(
  "C17_mCherry_FLAG_rep1",
  "C17_mCherry_NIPBL_rep1_rep2",
  "C17_mCherry_SMC1_rep2_rep3",
  "C17_mCherry_STAG1_rep1",
  "C17_mCherry_STAG2_rep1"
  )

y_axis_ranges <- list(
  "C17_mCherry_FLAG_rep1" = c(0, 7),
  "C17_mCherry_NIPBL_rep1_rep2" = c(0, 20),
  "C17_mCherry_SMC1_rep2_rep3" = c(0, 3),
  "C17_mCherry_STAG1_rep1" = c(0, 3),
  "C17_mCherry_STAG2_rep1" = c(0, 4)
)

mCherry_ASP <- plot_overview_grid(
  overview_exps = overview_exps,
  plots_by_expname = plots_by_expname,
  outF = outF,
  outfile_name = "250409_mCherry_ASP_overview_plot.pdf",
  title = "mCherry",
  y_axis_ranges = y_axis_ranges,
  remove_titles = TRUE
)

plot_overview_grid(
  overview_exps = overview_exps,
  plots_by_expname = plots_by_expname,
  outF = outF,
  outfile_name = "250409_mCherry_ASP_overview_plot_withTitle.pdf",
  title = "mCherry",
  y_axis_ranges = y_axis_ranges,
  remove_titles = FALSE
)

#SA2AID
overview_exps <- c(
  "C17_STAG2AID_C2_TACL-ON_IAA_FLAG_rep1_rep2", 
  "C17_STAG2AID_C2_TACL-ON_IAA_NIPBL_rep1",
  "C17_STAG2AID_C2_TACL-ON_IAA_SMC1_rep1", 
  "C17_STAG2AID_C2_TACL-ON_IAA_STAG1_rep1"
)

y_axis_ranges <- list(
  "C17_STAG2AID_C2_TACL-ON_IAA_FLAG_rep1_rep2" = c(0, 7),
  "C17_STAG2AID_C2_TACL-ON_IAA_NIPBL_rep1" = c(0, 20),
  "C17_STAG2AID_C2_TACL-ON_IAA_SMC1_rep1" = c(0, 3),
  "C17_STAG2AID_C2_TACL-ON_IAA_STAG1_rep1" = c(0, 3)
)

SA2_ASP <- plot_overview_grid(
  overview_exps = overview_exps,
  plots_by_expname = plots_by_expname,
  outF = outF,
  outfile_name = "250402_SA2AID_ASP_overview_plot.pdf",
  title = "SA2-AID IAA TACL-ON",
  y_axis_ranges = y_axis_ranges,
  remove_titles = TRUE
)

plot_overview_grid(
  overview_exps = overview_exps,
  plots_by_expname = plots_by_expname,
  outF = outF,
  outfile_name = "250402_SA2AID_ASP_overview_plot_withTitle.pdf",
  title = "SA2-AID IAA TACL-ON",
  y_axis_ranges = y_axis_ranges,
  remove_titles = FALSE
)





















#Aggregate 4C analysis ---------------------------------------------------------------









# 4C figures ---------------------------------------
#In 3C-based assays, ligation frequencies are typically highest near the viewpoint (<100 kb) and decrease as the distance from the viewpoint #increases. To minimize the high background ligation frequencies close to the viewpoint, only peaks located at least 100 kb away from the TetO #integration sites were included in the aggregate 4C analysis. These peaks were resized to 100 kb, divided into 2 kb bins, and the average normalized #4C signal was calculated for each bin.

#normalized_4C <- readRDS("/storage/shared/TACL/PK/4C/TACL_C17_degrons_normalized_4C.rds")

# Exclude the 17 and 19th TetO from the list as they have overlapping domains.
filtered_normalized_4C <- normalized_4C[-c(17, 19)]

# Combine the filtered list of GRanges objects into a single GRanges object
unlisted_normalized_4C <- do.call(c, filtered_normalized_4C)
start(unlisted_normalized_4C) <- unlisted_normalized_4C$pos
end(unlisted_normalized_4C) <- unlisted_normalized_4C$pos

#remove peaks within 100kb of the TetO
min_distance_Threshold <- 1e5
peaks_100kb <-peak_CTCF_C17_TMAU2[which(peak_CTCF_C17_TMAU2$distance > min_distance_Threshold)]

table(peaks_100kb$group, peaks_100kb$CTCF_strength_all, peaks_100kb$FIMO_convergent)





########Publication figure Fig 5 - 4C plots



getAvg4C_V2 <- function(region, ZoomWidth = 4e4, binSize = 1e3, unlisted_normalized_4C, colNames) {

  # # Process a single column:
  # tiles_out <- getAvg4C(region, unlisted_normalized_4C = unlisted_normalized_4C, colNames = "norm4C_on")

  # # Process multiple columns:
  # tiles_out <- getAvg4C(region, unlisted_normalized_4C = unlisted_normalized_4C, 
  #                       colNames = c("norm4C_on", "norm4C_off", "norm4C_STG2_on_IAA"))



  # Resize the region to the specified ZoomWidth (centered)
  ROI <- GenomicRanges::resize(region, width = ZoomWidth, fix = "center")
  # Tile the region into bins of size binSize
  ROI_tiles <- GenomicRanges::tile(x = ROI, width = binSize)
  tiles.gr <- unlist(ROI_tiles)
  
  # Find overlaps between the tiles and the 4C data
  hits <- findOverlaps(tiles.gr, unlisted_normalized_4C)
  hits_df <- as.data.frame(hits)
  
  # Loop over each specified column name
  for (colName in colNames) {
    # Extract the values for the current column from the 4C object
    values <- mcols(unlisted_normalized_4C)[[colName]]
    # Add these values to the hits data frame
    hits_df[[colName]] <- values[hits_df$subjectHits]
    
    # Compute the average signal per tile (query) for this column
    agg <- aggregate(hits_df[[colName]], 
                     by = list(query = hits_df$queryHits), 
                     FUN = mean)
    
    # Initialize the output metadata column for this signal (set to 0 by default)
    mcols(tiles.gr)[[colName]] <- numeric(length(tiles.gr))
    
    # Fill in the averaged values for each tile (using the tile indices from the aggregation)
    mcols(tiles.gr)[[colName]][agg$query] <- agg$x
  }
  
  return(tiles.gr)
}





plot_4C_allInOneGrid <- function(unlisted_normalized_4C, peaks,
                                 outF, outfile_name,
                                 ZoomWidth = 1e5, binSize = 2e3, 
                                 FlipRvStrand = FALSE,
                                 colNames = c("norm4C_on", "norm4C_off", "norm4C_mch",
                                              "norm4C_STAG2_on_IAA", "norm4C_WAPL_on_IAA",
                                              "norm4C_CTCF_on_IAA", "norm4C_PDS5A_on_IAA"),
                                 allowedPlotGroups = c("SA2collapsed_inside;convergent", 
                                                       "SA2collapsed_outside;convergent",
                                                       "SA2collapsed_inside;divergent", 
                                                       "SA2collapsed_outside;divergent"),
                                 group_col = "plotGroup_simple") {
  # Load necessary libraries
  library(GenomicRanges)
  library(ggplot2)
  library(tidyr)
  
  # Filter peaks based on the specified group column and non-missing CTCF_strength_all
  peaks_sub <- peaks[mcols(peaks)[, group_col] %in% allowedPlotGroups & 
                       !is.na(mcols(peaks)$CTCF_strength_all)]
  
  # Define the CTCF strength levels we want to plot (in a specific order)
  strength_levels <- c("Low", "Medium", "High")
  
  # Combined data frame for all conditions and plot groups
  agg_data_all <- data.frame()
  
  # Loop over each 4C condition (column)
  for (colName in colNames) {
    for (pg in allowedPlotGroups) {
      peaks_pg <- peaks_sub[mcols(peaks_sub)[, group_col] == pg]
      for (strength in strength_levels) {
        peaks_subset <- peaks_pg[mcols(peaks_pg)$CTCF_strength_all == strength]
        if (length(peaks_subset) == 0) next
        
        # Get binned coverage using your tiling function
        tiles <- getAvg4C_V2(region = peaks_subset,
                             ZoomWidth = ZoomWidth, binSize = binSize,
                             unlisted_normalized_4C = unlisted_normalized_4C,
                             colNames = colName)
        
        mat <- matrix(data = mcols(tiles)[, colName], 
                      nrow = length(peaks_subset), byrow = TRUE)
        
        if (FlipRvStrand) {
          idx_neg <- which(as.character(strand(peaks_subset)) == "-")
          if (length(idx_neg) > 0) {
            mat[idx_neg, ] <- t(apply(mat[idx_neg, , drop = FALSE], 1, rev))
          }
        }
        
        idx_no_signal <- which(rowSums(mat) == 0)
        if (length(idx_no_signal) > 0) {
          mat <- mat[-idx_no_signal, , drop = FALSE]
        }
        
        avg_signal <- colMeans(mat, na.rm = TRUE)
        n_bins <- length(avg_signal)
        x_positions <- seq(-ZoomWidth/2, ZoomWidth/2, length.out = n_bins)
        
        temp_df <- data.frame(
          Position   = x_positions,
          Signal     = avg_signal,
          Strength   = strength,
          PlotGroup  = pg,
          Condition  = colName,
          stringsAsFactors = FALSE
        )
        agg_data_all <- rbind(agg_data_all, temp_df)
      }
    }
  }
  
  if (nrow(agg_data_all) == 0) {
    message("No data aggregated.")
    return(NULL)
  }
  
  agg_data_all$Strength <- factor(agg_data_all$Strength, levels = strength_levels)
  agg_data_all$PlotGroup <- factor(agg_data_all$PlotGroup, levels = allowedPlotGroups)
  agg_data_all$Condition <- factor(agg_data_all$Condition, levels = colNames)
  
  # Define custom facet labels for Conditions
  condition_labels <- c(
    "norm4C_on"           = "TACL-ON",
    "norm4C_off"          = "TACL-OFF",
    "norm4C_mch"          = "Cherry",
    "norm4C_STAG2_on_IAA" = "STAG2 on IAA",
    "norm4C_WAPL_on_IAA"  = "WAPL on IAA",
    "norm4C_CTCF_on_IAA"  = "CTCF on IAA",
    "norm4C_PDS5A_on_IAA" = "PDS5A on IAA"
  )
  
  # Define custom labels for PlotGroup facets to match the ASP grid row order
  plotgroup_labels <- c(
    "SA2collapsed_inside;convergent" = "convergent\ninside",
    "SA2collapsed_outside;convergent" = "convergent\noutside",
    "SA2collapsed_inside;divergent"   = "divergent\ninside",
    "SA2collapsed_outside;divergent"   = "divergent\noutside"
  )
  
  # Build the combined plot with facet_grid: rows by PlotGroup, columns by Condition
  p <- ggplot(agg_data_all, aes(x = Position, y = Signal, color = Strength)) +
    geom_line(size = 1) +
    facet_grid(PlotGroup ~ Condition, 
               labeller = labeller(Condition = condition_labels, PlotGroup = plotgroup_labels)) +
    labs(title = "Aggregate 4C Signal",
         x = "Distance (bp)",
         y = "Average 4C Signal") +
    theme_minimal(base_size = 14) +
    theme(panel.grid.major = element_line(color = "gray80"),  # visible major grid lines
          panel.grid.minor = element_line(color = "gray90"),  # visible minor grid lines
          axis.line = element_line(color = "black"),          # x and y axis lines
          panel.border = element_rect(color = "black", fill = NA), # panel border
          plot.title = element_text(hjust = 0.5),
          aspect.ratio = 1) +
    scale_x_continuous(breaks = seq(-ZoomWidth/2, ZoomWidth/2, length.out = 5),
                       labels = round(seq(-ZoomWidth/2, ZoomWidth/2, length.out = 5)/1000, 0))
  
  # Set the output dimensions the same as your ASP overview grid:
  # Each Condition (column) gets 4 units of width plus an extra 1.5 units for the left margin.
  col_width <- 4
  width_total <- 1.5 + length(colNames) * col_width
  height_total <- 14  # Use the same overall height
  
  # Save the plot
  ggsave(file.path(outF, outfile_name),
         plot = p,
         width = width_total,
         height = height_total)
  
  return(p)
}



p_4C_grid <- plot_4C_allInOneGrid(
  unlisted_normalized_4C   = unlisted_normalized_4C,
  peaks                    = peaks_100kb,
  outF = outF,
  outfile_name = "agg4C_250402_100kbZoom_2kbBins_FlipTrue_12PM.pdf",
  ZoomWidth                = 1e5,
  binSize                  = 2e3,
  FlipRvStrand             = TRUE,
  colNames                 = c("norm4C_on", "norm4C_off", "norm4C_mch",
                               "norm4C_STAG2_on_IAA", "norm4C_WAPL_on_IAA",
                               "norm4C_CTCF_on_IAA", "norm4C_PDS5A_on_IAA"),
  allowedPlotGroups        = c("SA2collapsed_inside;convergent",
                               "SA2collapsed_outside;convergent",
                               "SA2collapsed_inside;divergent",
                               "SA2collapsed_outside;divergent")
)


