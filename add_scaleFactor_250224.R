
#R script to normalize ChIP-seq data for TACL manuscript revision
# ChIPdf_scalefactors_20250224.tsv is present in input folder Github page
# bigwigs can be downloaded from GEO. Modify the bigwig paths in the ChIPdf_scalefactors_20250224.tsv as needed.

# Load packages -----------------------------------------------------------

library("rtracklayer")
library("GenomicRanges")

# features ----------------------------------------------------------------

if (Sys.info()[[4]] == "srv-lnx-laat2") {
    prefix <- "/storage/shared/"
} else {
    prefix <- "/home/p.krijger_cbs-niob.local/mnt/laat2-storage/shared/"
}


C17TetO <- readRDS(file.path(prefix, "/TACL/PK/manuscript/revision/meta/final_C17locs_n27.rds"))
ChIPdf <- read.table("/storage/shared/TACL/PK/manuscript/revision/meta/ChIPdf_scalefactors_20250224.tsv", header = TRUE, sep = "\t")

# HMM domains
domain_HMM <- read.table("/storage/shared/TACL/PK/manuscript/revision/meta/TACL_domains_HMM.bed", header = FALSE, sep = "\t")
colnames(domain_HMM) <- c("chr", "start", "end", "name")
domain_HMM_GR <- GRanges(seqnames = domain_HMM$chr, ranges = IRanges(start = domain_HMM$start, end = domain_HMM$end))
# domain_HMM_GR$ID<-domain_HMM$name
domain_HMM_GR <- reduce(domain_HMM_GR) # remove redundant domains




# functions ----------------------------------------------------------------
source("./TACL_functions_240710.r")


# calculate scale fators ---------------------------------------------------


# Not normallized: H3K27me3, H3K4me1, H3K9me3


CTCF_norm_prots <- c("CTCF", "RAD21", "SMC1", "SMC3", "STAG1", "STAG2", "WAPL", "PDS5A")
FLAG_norm_prots <- c("FLAG", "MAU2", "NIPBL", "V5")

normProts <- c(CTCF_norm_prots, FLAG_norm_prots, "H3K27ac", "H3K4me3")

# add columns if not exist
if (!"meanBG" %in% colnames(ChIPdf)) {
    ChIPdf$meanBG <- NA
}

if (!"outerBG" %in% colnames(ChIPdf)) {
    ChIPdf$outerBG <- NA
}



for (a in 1:nrow(ChIPdf)) {
    if (anyNA(c(ChIPdf[a, "meanBG"], ChIPdf[a, "outerBG"]))) {
        message(ChIPdf$name[a])

        if (ChIPdf[a, "prot"] %in% normProts) {
            message("Normalizing")

            if (ChIPdf[a, "prot"] %in% CTCF_norm_prots) {
                peakFile <- "/TACL/PK/4DN_ChIP/outputs_VER9308/chip.peaks/eHAP1_C17_TMAU2_CTCF_rep1.idr.optimal_peak.regionPeak.bb"
                signalValueScore <- 35
            }

            if (ChIPdf[a, "prot"] %in% FLAG_norm_prots) {
                peakFile <- "/TACL/PK/4DN_ChIP/outputs_all/chip.peaks/TACL_RH_eHAP1_C17_T-hMAU2v2_FLAG_mrgd-VER6376_rep1-VER7254_rep2_input_VER6376_idr.optimal_peak.regionPeak.bb"
                signalValueScore <- 35
            }

            if (ChIPdf[a, "prot"] == "H3K27ac") {
                peakFile <- "/TACL/NKI/peaks/all_H3K27ac_peaks.canonical.replicated.no_blacklist.merged.bed"
                signalValueScore <- 0
            }

            if (ChIPdf[a, "prot"] == "H3K4me3") {
                peakFile <- "/TACL/NKI/peaks/all_H3K4me3_peaks.canonical.replicated.no_blacklist.merged.bed"
                signalValueScore <- 0
            }
        } else {
            message("Protein not in normProts. Skipping normalization.")
            next()
        }

        peakFile <- paste0(prefix, peakFile)

        if (!is.na(peakFile) & file.exists(peakFile)) {
            message(peakFile)
            peaks <- import(peakFile)
            peaks <- addTetOdistance(peaks, C17TetO)

            if (ChIPdf[a, "prot"] %in% c(CTCF_norm_prots, FLAG_norm_prots)) {
                peaks <- peaks[peaks$signalValue > 35, ]
                peaks <- filterOverlapPeaks(peaks = peaks)
            }
        } else {
            message("PeakFile not indicated or does not exist")
            next()
        }

        if (is.na(ChIPdf[a, "meanBG"])) {
            message("Calculating meanBG")
            meanBG <- getMeanBG_V3(region = peaks, bw = paste0(prefix, ChIPdf$bigwig[a]), signalValueScore = signalValueScore, distanceCutoff = 3e6, domain = domain_HMM_GR, binWidth = 10)
            message("meanBG:", meanBG)
            ChIPdf[a, "meanBG"] <- meanBG
        }

        if (is.na(ChIPdf[a, "outerBG"])) {
            message("Calculating outerBG")
            outerBG <- getOuterBinsBG_V3(region = peaks, bw = paste0(prefix, ChIPdf$bigwig[a]), signalValueScore = signalValueScore, distanceCutoff = 3e6, domain = domain_HMM_GR, zoomWidth = 5000, flankWidth = 1000)
            message("outerBG:", outerBG)
            ChIPdf[a, "outerBG"] <- outerBG
        }
    }
}


write.table(x = ChIPdf, file = paste0(prefix, "/TACL/PK/manuscript/revision/meta/ChIPdf_scalefactors_20250224.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)




