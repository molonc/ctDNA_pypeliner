require('tidyverse')
library('plyr')

read_results <- function(tsv_file) {
  results <- (read.csv(tsv_file, sep = "\t", header=TRUE, na.string = ".", quote = ""))
  if (nrow(results) > 0) {
    patient_id <- strsplit(tsv_file, "_")[[1]][1]
    file_name <- strsplit(tsv_file, "_")[[1]][2]
    file_name <- paste(strsplit(file_name, "-")[[1]][1], strsplit(file_name, "-")[[1]][2], sep="-")
    results$patient <- patient_id
    results$sample <- strsplit(file_name, "\\.")[[1]][1]
    results$chr <- as.character(results$chr)
    results$ref <- as.character(results$ref)
    results$alt <- as.character(results$alt)
    results$tn_ratio <- ifelse(results$N_vaf == 0, results$T_vaf * results$N_coverage, results$T_vaf/results$N_vaf)
    results$alt[results$alt == "TRUE"] <- "T"
    results$ref[results$ref == "TRUE"] <- "T"
    return(results)
  }
}

read_annotation <- function(tsv_file) {
  results <- (read.csv(tsv_file, sep = "\t", header=TRUE, na.string = ".", quote = ""))
  if (nrow(results) > 0) {
    patient_id <- strsplit(tsv_file, "_")[[1]][1]
    file_name <- strsplit(tsv_file, "_")[[1]][2]
    file_name <- paste(strsplit(file_name, "-")[[1]][1], strsplit(file_name, "-")[[1]][2], sep="-")
    results$patient <- patient_id
    results$sample <- strsplit(file_name, "\\.")[[1]][1]
    results$chr <- as.character(results$Chr)
    results$Ref <- as.character(results$Ref)
    results$Alt <- as.character(results$Alt)
    results$pos <- as.numeric(results$Start)
    results$pos <- ifelse(results$Alt == "-", results$pos - 1, results$pos)
    results$Alt[results$Alt == "TRUE"] <- "T"
    results$Ref[results$Ref == "TRUE"] <- "T"
    results <- subset(results, select = -c(Start, End, Chr))
    return(results)
  }
}

pointAlleleFrequency <- function(data){
  dataPlot <- ggplot(data = data, aes(x = N_vaf, y = T_vaf))
  dataPlot + geom_point(alpha = 0.8) + expand_limits(x = 0, y = 0)
}

plotByPosition <- function(data){
  ggplot(data, aes(pos, T_vaf, col = Gene.refGene))
  + geom_point() + facet_grid(~cytoBand, switch = 'x', scales = 'free_x')
  + theme(panel.spacing = unit(0.1, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_blank())
  + labs(x = 'cytoBand, pos')
}

plotByScore <- function(data){
  ggplot(plasma_snv_long, aes(value, T_vaf, col = as.factor(count)))
  + geom_point() + facet_grid(~tool, switch = 'x', scales = 'free')
  + theme(panel.spacing = unit(0.1, "lines"), strip.background = element_blank(), strip.placement = "outside")
  + labs(x = 'score')
}

histAlleleFrequency <- function(data){
  ggplot(data, aes(tn_ratio, fill = as.factor(count))) + geom_histogram() + scale_fill_brewer(palette = "Set3")
}

refAltDist <- function(data){
  ggplot(data = subset(data), aes(ref, fill=alt)) + geom_bar(position = "dodge")
}

filterSNVs <- function(data, match, tn_cutoff, match_tn_cutoff){
  data <- subset(data, !is.na(count))
  data <- merge(data, match, all.x = TRUE, by=c('chr','pos','patient','ref','alt'))
  colnames(data) <- sub("\\.x", "", colnames(data))
  data$match_tn_ratio <- data$tn_ratio.y
  data <- subset(data, (T_vaf >= 0.02 & tn_ratio >= tn_cutoff) | (!is.na(count.y) & match_tn_ratio >= match_tn_cutoff) | tn_ratio > 100)
  data <- data[, -grep("\\.y", colnames(data))]
  data <- subset(data, Func.refGene == 'exonic' & ExonicFunc.refGene != 'synonymous SNV')
  return(data)
}

filterIndels <- function(data){
  data <- subset(data, is.na(count))
  data <- subset(data, T_vaf >= 0.01 & tn_ratio >= 5)
  data <- subset(data, Func.refGene != 'intronic')
}

filterLOH <- function(data){
  data <- subset(data, N_vaf >= 0.25 & (tn_ratio >= 1.5 | tn_ratio <= 0.66))
  data <- subset(data, Func.refGene == 'exonic' & ExonicFunc.refGene != 'synonymous SNV')
}

annotations <- bind_rows(lapply(list.files(".", ".txt"), read_annotation))

plasma_validated <- subset(read.csv("plasma.csv", header=TRUE), vaf >= 1 | validation. == 'Confirmed')
tumour_validated <- subset(read.csv("tumour.csv", header=TRUE), vaf >= 1 & Sample.status != "WGA")
tumour_validated <- subset(tumour_validated, annotation_type == "conservative_inframe_insertion" | annotation_type == "frameshift_variant" | annotation_type == "missense_variant" | annotation_type == "stop_gained")
plasma_validated$validated <- TRUE
tumour_validated$validated <- TRUE

results <- bind_rows(lapply(list.files(".", ".tsv"), read_results))
results$type <- ifelse(startsWith(results$sample, "TNBC") | startsWith(results$sample, "VBA") | startsWith(results$sample, "PBC") | startsWith(results$sample, "BOB") | startsWith(results$sample, "TTR"), "plasma", "tumour")

results <- merge(results, annotations, by=c('chr', 'pos', 'patient', 'sample'))

plasma_results <- subset(results, type == "plasma")
tumour_results <- subset(results, type == "tumour")

validated_plasma_results <- merge(plasma_validated, plasma_results, all.x = TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
validated_tumour_results <- merge(tumour_validated, tumour_results, all.x = TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
plasma_tumour_snv <- merge(plasma_results, tumour_results, by=c('chr', 'pos', 'patient', 'ref', 'alt'))
plasma_results <- merge(plasma_results, plasma_validated, all.x=TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
tumour_results <- merge(tumour_results, tumour_validated, all.x=TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))

plasma_snv <- filterSNVs(plasma_results, tumour_results, 12, 20)
tumour_snv <- filterSNVs(tumour_results, plasma_results, 20, 12)

plasma_snv_long <- gather(plasma_snv, tool, value, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)
tumour_snv_long <- gather(tumour_snv, tool, value, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)

plasma_indel <- filterIndels(plasma_results)
tumour_indel <- filterIndels(tumour_results)

plasma_indel_long <- gather(plasma_indel, tool, value, VarScan, Strelka, factor_key = TRUE, na.rm = TRUE)
tumour_indel_long <- gather(tumour_indel, tool, value, VarScan, Strelka, factor_key = TRUE, na.rm = TRUE)

plasma_LOH <- filterLOH(plasma_results)
tumour_LOH <- filterLOH(tumour_results)

plasma_LOH_long <- gather(plasma_LOH, tool, vlaue, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)
tumour_LOH_long <- gather(tumour_LOH, tool, vlaue, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)
