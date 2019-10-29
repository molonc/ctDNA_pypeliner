setwd("~/shahlab/pye/projects/biof34/results/matched_results")
require('tidyverse')
library('plyr')

read_results <- function(tsv_file) {
  results <- (read.csv(tsv_file, sep = "\t", header=TRUE, na.string = ".", quote = ""))
  if (nrow(results) > 0) {
    patient_id <- strsplit(tsv_file, "_")[[1]][1]
    file_name <- strsplit(tsv_file, "_")[[1]][2]
    sample_name <- paste(sep="-", strsplit(file_name, "-")[[1]][1], strsplit(file_name, "-")[[1]][2])
    results$patient <- patient_id
    results$sample <- strsplit(sample_name, "\\.")[[1]][1]
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
    sample_name <- paste(sep="-", strsplit(file_name, "-")[[1]][1], strsplit(file_name, "-")[[1]][2])
    results$patient <- patient_id
    results$sample <- strsplit(sample_name, "\\.")[[1]][1]
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

classify_mutations <- function(data, match, vaf_cutoff, tn_cutoff, match_tn_cutoff){
  data <- merge(data, match, all.x = TRUE, by=c('chr','pos','patient','ref','alt'))
  colnames(data) <- sub("\\.x", "", colnames(data))
  data$match_tn_ratio <- data$tn_ratio.y
  data$classification <- ifelse((!is.na(data$count) & data$T_vaf >= vaf_cutoff & data$tn_ratio >= tn_cutoff), "high_freq_SNV", ifelse((!is.na(data$count.y) & data$match_tn_ratio >= match_tn_cutoff), "low_freq_in_tumour", NA))
  data$classification <- ifelse((is.na(data$count) & data$T_vaf >= vaf_cutoff & data$tn_ratio >= tn_cutoff), "indel", data$classification)
  data$classification <- ifelse((!is.na(data$count) & data$N_vaf >= 0.25 & (data$tn_ratio >= 1.5 | data$tn_ratio <= 0.66)), "LOH", data$classification)
  data <- data[, -grep("\\.y", colnames(data))]
  return(data)
}

## plot functions
# pointAlleleFrequency <- function(data){
#   dataPlot <- ggplot(data = data, aes(x = N_vaf, y = T_vaf))
#   dataPlot + geom_point(alpha = 0.8) + expand_limits(x = 0, y = 0)
# }
# 
# plotByPosition <- function(data){
#   ggplot(data, aes(pos, T_vaf, col = Gene.refGene)) + geom_point() + facet_grid(~cytoBand, switch = 'x', scales = 'free_x') + theme(panel.spacing = unit(0.1, "lines"), strip.background = element_blank(), strip.placement = "outside", axis.text.x = element_blank()) + labs(x = 'cytoBand, pos')
# }
# 
# plotByScore <- function(data){
#   ggplot(plasma_snv_long, aes(value, T_vaf, col = as.factor(count))) + geom_point() + facet_grid(~tool, switch = 'x', scales = 'free') + theme(panel.spacing = unit(0.1, "lines"), strip.background = element_blank(), strip.placement = "outside") + labs(x = 'score')
# }
# 
# histAlleleFrequency <- function(data){
#   ggplot(data, aes(tn_ratio, fill = as.factor(count))) + geom_histogram() + scale_fill_brewer(palette = "Set3")
# }
# 
# refAltDist <- function(data){
#   ggplot(data = data, aes(ref, fill=alt)) + geom_bar(position = "dodge")
# }

plot_patient_by_mutation <- function(data, type){
  data <- subset(data, classification != 'NA')
  events <- with(data, paste(Gene.refGene, AAChange.refGene))
  data$Gene.refGene <- factor(data$Gene.refGene, levels = rev(names(sort(table(sapply(strsplit(unique(events), split = " "), head, 1))))))
  data$AAChange.refGene <- sapply(strsplit(sapply(strsplit(as.character(data$AAChange.refGene), split = ","), head, 1), split = ":"), tail, 1)
  data$AAChange.refGene <- factor(data$AAChange.refGene, levels = rev(names(sort(table(data$AAChange.refGene)))))
  data$xaxis <- paste(data$Gene.refGene, data$AAChange.refGene)
  data$xaxis <- factor(data$xaxis, levels = unique(data$xaxis[order(data$Gene.refGene, data$AAChange.refGene)]))
  all_mutations <- ggplot(data, aes(xaxis, patient)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + aes(fill = ExonicFunc.refGene) + labs(x = "Gene + Amino Acid Change", fill = "Annotation Type", title = "Patient by Mutation Change", subtitle = "All Mutations") 
  snv_mutations <- ggplot(subset(data, classification == 'high_freq_SNV' | classification == 'low_freq_in_tumour'), aes(xaxis, patient)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + aes(fill = ExonicFunc.refGene) + labs(x = "Gene + Amino Acid Change", fill = "Annotation Type", title = "Patient by Mutation Change", subtitle = "SNV Mutations") 
  indel_mutations <- ggplot(subset(data, classification == 'indel'), aes(xaxis, patient)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + aes(fill = ExonicFunc.refGene) + labs(x = "Gene + Amino Acid Change", fill = "Annotation Type", title = "Patient by Mutation Change", subtitle = "Indel Mutations") 
  loh_mutations <- ggplot(subset(data, classification == 'LOH'), aes(xaxis, patient)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + aes(fill = ExonicFunc.refGene) + labs(x = "Gene + Amino Acid Change", fill = "Annotation Type", title = "Patient by Mutation Change", subtitle = "LOH Mutations") 
  pdf(paste0(type, "_patient_by_mutation.pdf"), width = 15, height = 12)
  print(all_mutations)
  print(snv_mutations)
  print(indel_mutations)
  # print(loh_mutations)
  dev.off()
}

plot_patient_plasma_timepoint <- function(){
  # plot patient plasma timepoint
  mapping <- read.csv("/Users/pye/pye_shah/ctDNA_pypeliner/helper_scripts/ctDNA.csv")
  mapping <- mapping[grepl('plasma', mapping$Type),]
  mapping$Type <- factor(mapping$Type)
  mapping$sample <- with(mapping, gsub(" ", "-", paste(Aliquot.ID, Alias)))
  mapping$Type <- gsub('plasma ', '', mapping$Type)
  mapping$Type <- gsub(' plasma', '', mapping$Type)
  # mapping$Type <- gsub("19", "18", mapping$Type)
  # mapping$Type <- gsub("16", "15", mapping$Type)
  # mapping$Type <- gsub("earliest ", "", mapping$Type)
  # mapping$Type <- gsub("month repeat", "month", mapping$Type)
  # mapping <- subset(mapping, Type != 'WGA baseline')
  mapping$status <- ifelse(mapping$sample %in% unique(subset(plasma, classification != 'NA')$sample), 'ctDNA +', 'ctDNA -')
  mapping$Type <- factor(mapping$Type, levels = c('earliest baseline', 'baseline', 'baseline repeat', 'WGA baseline', '1 month', '3 month', '6 month', '6 month repeat', '9 month', '12 month', '12 month repeat', '15 month', '16 month', '18 month', '18 month repeat', '19 month', '21 month', '24 month', '30 month', '36 month'))
  patient_mapping <- read.csv("patient_mapping.csv", header = TRUE, na.string = c("N/A", ''))
  patient_mapping <- gather(patient_mapping, Type, Date, baseline:X36.month, factor_key = TRUE)
  patient_mapping$Type <- gsub("X", "", gsub("\\.", " ", patient_mapping$Type))
  masterlist <- merge(patient_mapping, mapping, all.x = TRUE, all.y = TRUE, by = c('PBC.ID', 'Type'))
  masterlist$status <- ifelse(is.na(masterlist$status) & !is.na(masterlist$Date), "collected, not sequenced", masterlist$status)
  masterlist$Type <- factor(masterlist$Type, levels = levels(mapping$Type))
  masterlist$PBC.ID <- as.character(masterlist$PBC.ID)
  pdf("patient_plasma_timepoint.pdf", width = 15, height = 12)
  print(ggplot(masterlist, aes(Type, PBC.ID, fill = status)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("#F7BB8C", "#5DDDB9", "#FF8F8F"), na.value="#D6E3F8"))
  dev.off()
}

plot_concordance <- function(plasma, tumour){
  # plot concordance/discordance between plasma and tumour
  merge_results <- merge(subset(plasma, classification != 'NA'), subset(tumour, classification != 'NA'), all = TRUE, by = c('patient', 'ref', 'alt', 'chr', 'pos', 'cytoBand'))
  merge_results$presence <- ifelse(!is.na(merge_results$sample.x) & !is.na(merge_results$sample.y), "Both", ifelse(!is.na(merge_results$sample.x), "Plasma Only", "Tumour Only"))
  pdf("plasma_tumour_concordance.pdf", width = 15, height = 12)
  print(ggplot(merge_results, aes(patient, fill = presence)) + geom_bar(data = subset(merge_results, presence == "Both"), aes(y = stat(count))) + geom_bar(data = subset(merge_results, presence == "Plasma Only" | presence == "Tumour Only"), aes(y = -stat(count)), position = "stack") + coord_flip() + scale_y_continuous(breaks = seq(-40, 50, 10), labels = paste0(as.character(c(4:0, 1:5)), "0")))
  dev.off()
}

annotations <- bind_rows(lapply(list.files(".", ".txt"), read_annotation))

# validated data set
plasma_validated <- subset(read.csv("plasma.csv", header=TRUE), validation. == 'Confirmed')
plasma_validated$patient <- plasma_validated$PBC.ID
plasma_validated <- subset(plasma_validated, select = c(chr, pos, sample, patient, ref, alt, validation.))

# read results
results <- bind_rows(lapply(list.files(".", ".tsv"), read_results))
results$type <- ifelse(startsWith(results$sample, "TNBC") | startsWith(results$sample, "VBA") | startsWith(results$sample, "PBC") | startsWith(results$sample, "BOB") | startsWith(results$sample, "TTR"), "plasma", "tumour")

results <- merge(results, annotations, by=c('chr', 'pos', 'patient', 'sample'))

plasma_results <- subset(results, type == "plasma")
tumour_results <- subset(results, type == "tumour")
plasma_tumour <- merge(plasma_results, tumour_results, by=c('chr', 'pos', 'patient', 'ref', 'alt'))

# compare against validated
plasma_results <- merge(plasma_results, plasma_validated, all.x=TRUE, by=c('chr', 'pos', 'sample', 'patient', 'ref', 'alt'))

# classify mutation type
plasma <- classify_mutations(plasma_results, tumour_results, 0.02, 10, 15)
tumour <- classify_mutations(tumour_results, plasma_results, 0.02, 15, 10)

# plot results
plot_patient_by_mutation(plasma, "plasma")
plot_patient_by_mutation(tumour, "tumour")

plot_patient_plasma_timepoint()

plot_concordance(plasma, tumour)

# separate by classification
plasma_snv <- subset(plasma, classification == 'high_freq_SNV' | classification == 'low_freq_in_tumour')
plasma_indel <- subset(plasma, classification == 'indel')
plasma_LOH <- subset(plasma, classification == 'LOH')

tumour_snv <- subset(tumour, classification == 'high_freq_SNV' | classification == 'low_freq_in_tumour')
tumour_indel <- subset(tumour, classification == 'indel')
tumour_LOH <- subset(tumour, classification == 'LOH')

plasma_snv <- subset(plasma_snv, select = -c(T_ref, T_alt, N_ref, N_alt))
tumour_snv <- subset(tumour_snv, select = -c(T_ref, T_alt, N_ref, N_alt))

plasma_indel <- subset(plasma_indel, select = -c(count, deepSNV, LoLoPicker, MutationSeq, N_A, N_C, N_G, N_T, N_N, T_A, T_C, T_G, T_T, T_N))
tumour_indel <- subset(tumour_indel, select = -c(count, deepSNV, LoLoPicker, MutationSeq, N_A, N_C, N_G, N_T, N_N, T_A, T_C, T_G, T_T, T_N))

plasma_LOH <- subset(plasma_LOH, select = -c(T_ref, T_alt, N_ref, N_alt))
tumour_LOH <- subset(tumour_LOH, select = -c(T_ref, T_alt, N_ref, N_alt))

# separate by tool
plasma_long <- gather(plasma, tool, value, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)
tumour_long <- gather(tumour, tool, value, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)
plasma_snv_long <- gather(plasma_snv, tool, value, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)
tumour_snv_long <- gather(tumour_snv, tool, value, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)
plasma_indel_long <- gather(plasma_indel, tool, value, VarScan, Strelka, factor_key = TRUE, na.rm = TRUE)
tumour_indel_long <- gather(tumour_indel, tool, value, VarScan, Strelka, factor_key = TRUE, na.rm = TRUE)
plasma_LOH_long <- gather(plasma_LOH, tool, value, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)
tumour_LOH_long <- gather(tumour_LOH, tool, value, deepSNV, LoLoPicker, VarScan, MutationSeq, Strelka, factor_key = TRUE, na.rm = TRUE)

# write to CSV files
write.csv(plasma_indel, "plasma_indel.csv", row.names = FALSE)
write.csv(tumour_indel, "tumour_indel.csv", row.names = FALSE)
write.csv(plasma_LOH, "plasma_LOH.csv", row.names = FALSE)
write.csv(tumour_LOH, "tumour_LOH.csv", row.names = FALSE)
write.csv(plasma_snv, "plasma_snv.csv", row.names = FALSE)
write.csv(tumour_snv, "tumour_snv.csv", row.names = FALSE)
write.csv(plasma, "plasma_unfiltered.csv", row.names = FALSE)
write.csv(tumour, "tumour_unfiltered.csv", row.names = FALSE)
