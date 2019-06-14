#!/usr/bin/env Rscript
require("deepSNV")
require("optparse")

option_list <- list(make_option(c("-n", "--normal"), type="character", default=NULL, help="normal_bam", metavar="character"),
                    make_option(c("-t", "--tumour"), type="character", default=NULL, help="tumour_bam", metavar="character"),
                    make_option(c("-b", "--bed"), type="character", default="~/BIOF-34/ctDNA_pypeliner/beds/CG001v4.0.bed", help="bed_file", metavar="character"),
                    make_option(c("-q", "--quality"), type="integer", default=25, help="phred_quality", metavar="character"),
                    make_option(c("-o", "--out"), type="character", default="deepSNV_out.tsv", help="output_file", metavar="character")
                    )

opt_parser <- OptionParser(option_list=option_list)

opt <- parse_args(opt_parser)

bed <- read.table(opt$bed, header=FALSE, sep="\t", quote="", stringsAsFactors=FALSE)
regions <- data.frame(chr=bed[,"V1"], start=bed[,"V2"], stop=bed[,"V3"])

analysis <- deepSNV(test=opt$tumour, control=opt$normal, q=opt$quality, regions=regions)
summary <- summary(analysis)

write.table(summary, opt$out, sep="\t", row.names=FALSE, quote=FALSE)
