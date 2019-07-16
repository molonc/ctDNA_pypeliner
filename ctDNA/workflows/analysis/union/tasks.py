from __future__ import division
import csv
import vcf
from pypeliner.commandline import execute

def bam_readcount(config, bam_type, bam, result, tmp_file):
    execute(
        'bam-readcount',
        '-d',
        20000,
        '-w',
        0,
        '--min-mapping-quality',
        50,
        '--min-base-quality',
        30,
        '-f',
        config['reference_genome'],
        bam,
        result['chr'] + ":" + result['pos'] + "-" + result['pos'],
        '>',
        tmp_file,
        )

    with open(tmp_file, 'rb') as temp:
        for row in temp:
            read_counts = row.split()
            result[bam_type + '_coverage'] = float(read_counts[3])
            for nucleotide_count in read_counts[5:10]:
                nucleotide_count_split = nucleotide_count.split(":")
                result[bam_type + '_' + nucleotide_count_split[0]] = float(nucleotide_count_split[1])

            result[bam_type + '_alf'] = result[bam_type + '_' + result['alt']]/result[bam_type + '_coverage']

def DeepSNV_process(tool, input_file, results):
    with open(input_file, 'rb') as inputs:
        reader = csv.DictReader(inputs, delimiter='\t')
        for row in reader:
            key = row['chr'] + ':' + row['pos']
            try:
                taf = round(float(row['freq.var']), 4)
            except ZeroDivisionError:
                taf = 0

            if taf < 0.001 or row['var'] == '-' or float(row['p.val']) > 0.0005:
                continue

            if results.get(key, False) and not results[key].get(tool, False):
                results[key]['count'] += 1
                results[key][tool] = row['p.val']

            else:
                results[key] = {
                    'chr': row['chr'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['var'],
                    'count': 1,
                    tool: row['p.val']
                    }

def LoLoPicker_process(tool, input_file, results):
    with open(input_file, 'rb') as inputs:
        reader = csv.DictReader(inputs, delimiter='\t')
        for row in reader:
            key = row['#chr'] + ':' + row['pos']
            try:
                taf = round(float(row['tumor_alf']), 4)
            except ZeroDivisionError:
                taf = 0

            if taf < 0.001:
                continue

            if results.get(key, False) and not results[key].get(tool, False):
                results[key]['count'] += 1
                results[key][tool] = row['p_value']

            else:
                results[key] = {
                    'chr': row['#chr'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['alt'],
                    'count': 1,
                    tool: row['p_value']
                    }

def VarScan_process(tool, input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        freq = row.genotype('TUMOR')['FREQ']
        try:
            taf = round(float(freq[:-1]) / 100, 4)
        except ZeroDivisionError:
            taf = 0

        if taf < 0.001 or float(row.INFO['SPV']) > 0.0005:
            continue

        if results.get(key, False) and not results[key].get(tool, False):
            results[key]['count'] += 1
            results[key][tool] = row.INFO['SPV']

        else:
            for alt in row.ALT:
                results[key] = {
                    'chr': row.CHROM,
                    'pos': str(row.POS),
                    'ref': str(row.REF),
                    'alt': str(alt),
                    'count': 1,
                    tool: str(row.INFO['SPV'])
                    }

def museq_process(tool, input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        try:
            taf = round(float(row.INFO['TA']) / (float(row.INFO['TA']) + float(row.INFO['TR'])), 4)
        except ZeroDivisionError:
            taf = 0

        if taf < 0.001 or float(row.INFO['PR']) < 0.65:
            continue

        if results.get(key, False) and not results[key].get(tool, False):
            results[key]['count'] += 1
            results[key][tool] = str(row.INFO['PR'])

        else:
            for alt in row.ALT:
                results[key] = {
                    'chr': row.CHROM,
                    'pos': str(row.POS),
                    'ref': str(row.REF),
                    'alt': str(alt),
                    'count': 1,
                    tool: str(row.INFO['PR'])
                    }

def strelka_process(tool, input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        try:
            taf = round(float(row.genotype('TUMOR')[str(row.ALT[0]) + 'U'][0])/float(row.genotype('TUMOR')['DP']), 4)
        except ZeroDivisionError:
            taf = 0

        if taf < 0.001 or int(row.INFO['QSS']) < 200:
            continue

        if results.get(key, False) and not results[key].get(tool, False):
            results[key]['count'] += 1
            results[key][tool] = str(row.INFO['QSS'])

        else:
            for alt in row.ALT:
                results[key] = {
                    'chr': row.CHROM,
                    'pos': str(row.POS),
                    'ref': str(row.REF),
                    'alt': str(alt),
                    'count': 1,
                    tool: str(row.INFO['QSS'])
                    }