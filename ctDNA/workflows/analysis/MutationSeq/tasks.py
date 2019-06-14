import os
import csv
import pysam
from pypeliner.commandline import execute
from ctDNA.utils import vcfutils

def generate_intervals(ref, size=1000000):
    chromosomes = list(map(lambda x: 'chr' + str(x), range(1, 23) + ['X']))
    fasta = pysam.FastaFile(ref)
    lengths = fasta.lengths
    print lengths
    names = fasta.references
    print names

    intervals = []

    for name, length in zip(names, lengths):
        if name not in chromosomes:
            continue
        for i in range(int((length / size) + 1)):
            start = str(int(i * size))
            end = str(int((i + 1) * size))
            intervals.append(name + "_" + start + "_" + end)

    print intervals
    return intervals

def run_museq(config, normal_bam, tumour_bam, interval, output_file, log_file):
    interval = interval.split('_')
    interval = interval[0] + ':' + interval[1] + '-' + interval[2]

    execute(
        'museq',
        'normal:' + normal_bam,
        'tumour:' + tumour_bam,
        'reference:' + config['reference_genome'],
        '--deep',
        '--manifest',
        config['bed_file'],
        '--interval',
        interval,
        '--out',
        output_file,
        '--log',
        log_file
        )

def merge_vcfs(inputs, outfile, tempdir):
    helpers.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, outfile)
