import os
import csv
from pypeliner.commandline import execute
from ctDNA.utils import vcfutils

def generate_intervals(bed_file):
	intervals = []
	with open(bed_file, 'r') as bed:
		bed_reader = csv.reader(bed, delimiter='\t')
		for row in bed_reader:
			intervals.append(row[0] + ':' + row[1] + '-' + row[2])

	return intervals

def run_museq(config, normal_bam, tumour_bam, interval, output_file, log_file):
	execute(
		'museq',
		'normal:' + normal_bam,
		'tumour:' + tumour_bam,
		'reference:' + config['reference_genome'],
		'--interval',
		interval,
		'--out',
		output_file,
		'--log',
		log_file
		)

def merge_vcfs(inputs, outfile):
	vcfutils.concatenate_vcf(inputs, outfile)
