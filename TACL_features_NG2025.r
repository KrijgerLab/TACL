

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


# Features ----------------------------------------------------------------

prefix <- ""
folder_4C <- "./normalized_4c"   #Github page
folder_peaks <- "./peaks"  #Github page

# load TetO site locations
    # C17TetO <- read.table(file.path(prefix,'/TACL/PK/C17/final_C17locs_n27.bed'))
    # colnames(C17TetO) <- c('seqnames','start','end','score4C','judgement','TetO')
    # C17TetO <- GRanges(C17TetO)
    # saveRDS(C17TetO, file = '/storage/shared/./input/final_C17locs_n27.rds')
C17TetO <- readRDS(file.path(prefix, "./input/final_C17locs_n27.rds"))
protCols <- read.table(file.path(prefix, "./input/protCol_manuscript.tsv"), header = TRUE, comment.char = "")

dCols <- protCols[protCols$prot %in% c("FLAG", "MAU2", "SMC1", "RAD21"), ]
dCols$protein <- paste0("d", dCols$protein)
rbind(protCols, dCols)



# ChIP info
ChIPdf <- read.table(file.path(prefix, "./input/ChIPdf_scalefactors_20250224.tsv"), header = TRUE, sep = "\t")

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
TetO_enrichment <- read.table(file.path(prefix, "./input/TetO_enrichment_250224.tsv"), header = TRUE, sep = "\t")
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
domain_HMM <- read.table(file.path(prefix, "./input/TACL_domains_HMM.bed"), header = FALSE, sep = "\t")
colnames(domain_HMM) <- c("chr", "start", "end", "name")
domain_HMM_GR <- GRanges(seqnames = domain_HMM$chr, ranges = IRanges(start = domain_HMM$start, end = domain_HMM$end))
# domain_HMM_GR$ID<-domain_HMM$name  #note mikhail added the wrong name for 2 and 25
domain_HMM_GR <- reduce(domain_HMM_GR) # remove redundant domains


# CTCF
peak_CTCF_C17_TMAU2 <- readRDS(file.path(prefix, "peaks/C17_TACL_CTCF_hg38_peaks_withMotif_TetO.rds"))
peak_CTCF_C17_TMAU2$HMMdomain <- countOverlaps(peak_CTCF_C17_TMAU2, domain_HMM_GR)

peak_CTCF_C17_Tmcherry <- readRDS(file.path(prefix, "peaks/C17_mCherry_CTCF_hg38_peaks_withMotif_TetO.rds"))


CTCF_reduced <- reduce(c(peak_CTCF_C17_TMAU2, peak_CTCF_C17_Tmcherry))

CTCF_reduced_sV35 <- reduce(c(peak_CTCF_C17_TMAU2[peak_CTCF_C17_TMAU2$signalValue > 35], peak_CTCF_C17_Tmcherry[peak_CTCF_C17_Tmcherry$signalValue > 35]))



peak_CTCF_C17_TMAU2_convergent <- peak_CTCF_C17_TMAU2[which(peak_CTCF_C17_TMAU2$FIMO_convergent)] # 14303
peak_CTCF_C17_TMAU2_divergent <- peak_CTCF_C17_TMAU2[which(peak_CTCF_C17_TMAU2$FIMO_convergent == FALSE)] # 14210



# C17 Flag peaks in TACL-ON
flag_opt <- import(file.path(prefix, "peaks/TACL_RH_eHAP1_C17_T-hMAU2v2_FLAG_mrgd-VER6376_rep1-VER7254_rep2_input_VER6376_idr.optimal_peak.regionPeak.bb")) # 11757
flag_opt_mrgdbyhighestPeak <- filterOverlapPeaks(peaks = flag_opt) # 10442
flag_opt_mrgdbyhighestPeak$ID <- 1:length(flag_opt_mrgdbyhighestPeak)
flag_opt_mrgdbyhighestPeak <- addTetOdistance(flag_opt_mrgdbyhighestPeak, C17TetO)
flag_opt_mrgdbyhighestPeak$HMMdomain <- countOverlaps(flag_opt_mrgdbyhighestPeak, domain_HMM_GR)

FLAG_filtered <- flag_opt_mrgdbyhighestPeak[flag_opt_mrgdbyhighestPeak$signalValue > 35] # 5725 peaks
table(FLAG_filtered$distance < 3e6) # F 4795; T 930

# dFLAG peaks
dFLAG <- readRDS(file.path(prefix, "./peaks/dFLAG_peaks.rds"))

# NKI peaks
NKI_H3K27ac_peaks <- import(file.path(prefix, "/peaks/all_H3K27ac_peaks.canonical.replicated.no_blacklist.merged.bed"))
NKI_H3K4me3_peaks <- import(file.path(prefix, "/peaks/all_H3K4me3_peaks.canonical.replicated.no_blacklist.merged.bed"))

# H3K4me1 (>=1)
H3K4me1_opt_All <- import(paste0(prefix, "/peaks/Encode_H3K4me1/HAP1_H3K4me1_ENCODE_hg38_pseudorep_ENCFF963OJW.bb"))
H3K4me1_opt_All <- H3K4me1_opt_All[H3K4me1_opt_All$signalValue >= 1]
H3K4me1_opt_reduced <- GenomicRanges::reduce(H3K4me1_opt_All) # 116897



# STAG2 lost domains
STAG2_lost_domains <- read.table("./input/STAG2_lost_domains.txt", header = TRUE, sep = "\t")
STAG2_lost_domains_GR <- GRanges(seqnames = STAG2_lost_domains$chrom, ranges = IRanges(start = STAG2_lost_domains$tacl_domain_stag2_iaa_start, end = STAG2_lost_domains$tacl_domain_stag2_iaa_end))
STAG2_lost_domains_GR$teto_id <- STAG2_lost_domains$teto_id

STAG2_lost_domains_WT_GR <- GRanges(seqnames = STAG2_lost_domains$chrom, ranges = IRanges(start = STAG2_lost_domains$tacl_domain_stag2_wt_start, end = STAG2_lost_domains$tacl_domain_stag2_wt_end))
STAG2_lost_domains_WT_GR$teto_id <- STAG2_lost_domains$teto_id


# degron domains
degron_domains <- read.table("./input/degrons_domains_HMM.txt", header = TRUE, sep = "\t")
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
# normalized_4C_rep2 <- readRDS(file.path(folder_4C, "TACL_rep2_normalized_4C.rds"))
# normalized_4C <- readRDS(file.path(folder_4C,"TACL_C17_degrons_including_mch_VER10253_normalized_4C.rds"))
# normalized_4C_RAD21 <- readRDS(file.path(folder_4C,"TACL_C17_RAD21AID_normalized_4C_240709.rds"))




