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
    results$ref <- as.character(results$Ref)
    results$alt <- as.character(results$Alt)
    results$pos <- as.numeric(results$Start)
    results$pos <- ifelse(results$alt == "-", results$pos - 1, results$pos)
    results$alt[results$alt == "TRUE"] <- "T"
    results$ref[results$ref == "TRUE"] <- "T"
    results <- subset(results, select = -c(Start, End, Chr, Ref, Alt))
    return(results)
  }
}

pointAlleleFrequency <- function(data){
  dataPlot <- ggplot(data = data, aes(x = N_vaf, y = T_vaf))
  dataPlot + geom_point(alpha = 0.8) + expand_limits(x = 0, y = 0)
}

pointValidated <- function(data){
  tn_ratio_mean <- mean(data$tn_ratio)
  std_dev <- sd(data$tn_ratio)
  data <- subset(data, tn_ratio > 1.5)
  dataPlot <- ggplot(data = data, aes(x = N_vaf, y = T_vaf, col = as.factor(count)))
  unlabeled <- dataPlot + geom_point(alpha = 0.8) + geom_abline(slope = 1) + scale_color_brewer(palette = "Set3")
  unlabeled + geom_text(aes(label=ifelse(validation. == "Confirmed", paste(chr, ":", pos), '')), hjust=0, vjust=0)  
}

histAlleleFrequency <- function(data){
  ggplot(data, aes(tn_ratio, fill = as.factor(count))) + geom_histogram() + scale_fill_brewer(palette = "Set3")
}

refAltDist <- function(data){
  ggplot(data = subset(data), aes(ref, fill=alt)) + geom_bar(position = "dodge")
}

annotations <- bind_rows(lapply(list.files(".", ".txt"), read_annotation))

SNVs <- bind_rows(lapply(list.files(".", ".snv.tsv"), read_results))
SNVs$type <- ifelse(startsWith(SNVs$sample, "TNBC") | startsWith(SNVs$sample, "TTR") | startsWith(SNVs$sample, "PBC") | startsWith(SNVs$sample, "BOB"), "plasma", "tumour")
plasma_snv <- subset(SNVs, type == "plasma")
tumour_snv <- subset(SNVs, type == "tumour")
# tumour_snv <- subset(tumour_snv, !((ref == "C" & alt == "T") | (ref == "G" & alt == "A")) | (!is.na(Strelka) & !is.na(MutationSeq)))

plasma_validated <- read.csv("plasma.csv", header=TRUE)
tumour_validated <- read.csv("tumour.csv", header=TRUE)
plasma_validated$validated <- TRUE
tumour_validated$validated <- TRUE
plasma_validated <- subset(plasma_validated, mutation_type == 'snv' & (vaf >= 1 | validation. == 'Confirmed'))
tumour_validated <- subset(tumour_validated, mutation_type == 'snv' & vaf >= 1 & Sample.status != "WGA")
tumour_validated <- subset(tumour_validated, annotation_type == "missense_variant" | annotation_type == "stop_gained")

validated_plasma_snv <- merge(plasma_validated, plasma_snv, all.x = TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
validated_tumour_snv <- merge(tumour_validated, tumour_snv, all.x=TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
plasma_snv <- merge(plasma_snv, plasma_validated, all.x=TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
tumour_snv <- merge(tumour_snv, tumour_validated, all.x=TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
plasma_tumour_snv <- merge(plasma_snv, tumour_snv, by=c('chr', 'pos', 'patient', 'ref', 'alt'))

SNVs_long <- gather(SNVs, tool, value, deepSNV:Strelka, factor_key=TRUE)
SNVs_long <- subset(SNVs_long, !is.na(value))
plasma_snv_long <- gather(plasma_snv, tool, value, deepSNV:Strelka, factor_key = TRUE)
plasma_snv_long <- subset(plasma_snv_long, !is.na(value))
tumour_snv_long <- gather(tumour_snv, tool, value, deepSNV:Strelka, factor_key = TRUE)
tumour_snv_long <- subset(tumour_snv_long, !is.na(value))

plasma_snvs <- merge(plasma_snv, tumour_snv, all.x = TRUE, by=c('chr','pos','patient','ref','alt'))
colnames(plasma_snvs) <- sub("\\.x", "", colnames(plasma_snvs))
plasma_snvs$tumour_tn_ratio <- plasma_snvs$tn_ratio.y
plasma_snvs <- subset(plasma_snvs, (T_vaf >= 0.02 & tn_ratio >= 12) | (!is.na(count.y) & tumour_tn_ratio >= 20))
plasma_snvs <- plasma_snvs[, -grep("\\.y", colnames(plasma_snvs))]
plasma_snvs <- merge(plasma_snvs, annotations, by=c('chr', 'pos', 'patient', 'ref', 'alt', 'sample'))
tumour_snvs <- merge(tumour_snv, plasma_snv, all.x = TRUE, by=c('chr','pos','patient','ref','alt'))
colnames(tumour_snvs) <- sub("\\.x", "", colnames(tumour_snvs))
tumour_snvs$plasma_tn_ratio <- tumour_snvs$tn_ratio.y
tumour_snvs <- subset(tumour_snvs, (T_vaf >= 0.02 & tn_ratio >= 20) | (!is.na(count.y) & plasma_tn_ratio >= 12))
tumour_snvs <- tumour_snvs[, -grep("\\.y", colnames(tumour_snvs))]
tumour_snvs <- merge(tumour_snvs, annotations, by=c('chr', 'pos', 'patient', 'ref', 'alt', 'sample'))

plasma_snvs <- subset(plasma_snvs, Func.refGene == 'exonic' & ExonicFunc.refGene != 'synonymous SNV')
tumour_snvs <- subset(tumour_snvs, Func.refGene == 'exonic' & ExonicFunc.refGene != 'synonymous SNV')

plasma_snvs_long <- gather(plasma_snvs, tool, value, deepSNV:Strelka, factor_key = TRUE)
plasma_snvs_long <- subset(plasma_snvs_long, !is.na(value))
tumour_snvs_long <- gather(tumour_snvs, tool, value, deepSNV:Strelka, factor_key = TRUE)
tumour_snvs_long <- subset(tumour_snvs_long, !is.na(value))

indels <- bind_rows(lapply(list.files(".", ".indel.tsv"), read_results))
indels$type <- ifelse(startsWith(indels$sample, "TNBC") | startsWith(indels$sample, "VBA") | startsWith(indels$sample, "PBC") | startsWith(indels$sample, "BOB") | startsWith(indels$sample, "TTR"), "plasma", "tumour")
plasma_indel <- subset(indels, type == "plasma")
tumour_indel <- subset(indels, type == "tumour")

plasma_validated <- read.csv("plasma.csv", header=TRUE)
tumour_validated <- read.csv("tumour.csv", header=TRUE)
plasma_validated$validated <- TRUE
tumour_validated$validated <- TRUE
plasma_validated <- subset(plasma_validated, mutation_type == 'indel' & (vaf >= 1 | validation. == 'Confirmed'))
tumour_validated <- subset(tumour_validated, mutation_type == 'indel' & vaf >= 1 & Sample.status != "WGA")

plasma_indel <- merge(plasma_indel, plasma_validated, all.x=TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
tumour_indel <- merge(tumour_indel, tumour_validated, all.x=TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
plasma_tumour_indel <- merge(tumour_indel, plasma_indel, by=c('chr', 'pos', 'patient', 'ref', 'alt'))
validated_plasma_indel <- merge(plasma_validated, plasma_indel, all.x = TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))
validated_tumour_indel <- merge(tumour_validated, tumour_indel, all.x=TRUE, by=c('chr', 'pos', 'sample', 'ref', 'alt'))

indels_long <- gather(indels, tool, value, VarScan:Strelka, factor_key=TRUE)
indels_long <- subset(indels_long, !is.na(value))
plasma_indel_long <- gather(plasma_indel, tool, value, VarScan:Strelka, factor_key = TRUE)
plasma_indel_long <- subset(plasma_indel_long, !is.na(value))
tumour_indel_long <- gather(tumour_indel, tool, value, VarScan:Strelka, factor_key = TRUE)
tumour_indel_long <- subset(tumour_indel_long, !is.na(value))

plasma_indels <- subset(plasma_indel, T_vaf >= 0.01 & tn_ratio >= 5)
plasma_indels <- merge(plasma_indels, annotations, by=c('chr', 'patient', 'sample', 'pos'))
tumour_indels <- subset(tumour_indel, T_vaf >= 0.01 & tn_ratio >= 5)
tumour_indels <- merge(tumour_indels, annotations, by=c('chr', 'patient', 'sample', 'pos'))

plasma_indels_long <- gather(plasma_indels, tool, value, VarScan:Strelka, factor_key = TRUE)
plasma_indels_long <- subset(plasma_indels_long, !is.na(value))
tumour_indels_long <- gather(tumour_indels, tool, value, VarScan:Strelka, factor_key = TRUE)
tumour_indels_long <- subset(tumour_indels_long, !is.na(value))