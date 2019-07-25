from __future__ import division
import csv
import vcf
from pypeliner.commandline import execute

T_VAF_CUTOFF = 0.001
P_VALUE_CUTOFF = 0.0005
MUSEQ_CUTOFF = 0.65
STRELKA_CUTOFF = 200

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

            result[bam_type + '_vaf'] = result[bam_type + '_' + result['alt']]/result[bam_type + '_coverage']

def DeepSNV_process(tool, input_file, results):
    with open(input_file, 'rb') as inputs:
        reader = csv.DictReader(inputs, delimiter='\t')
        for row in reader:
            key = row['chr'] + ':' + row['pos']
            t_vaf = round(float(row['freq.var']), 4)

            if t_vaf < 0.001 or row['var'] == '-' or float(row['p.val']) > P_VALUE_CUTOFF:
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
            t_vaf = round(float(row['tumor_alf']), 4)

            if t_vaf < T_VAF_CUTOFF or float(row['corrected_p']) > P_VALUE_CUTOFF:
                continue

            if results.get(key, False) and not results[key].get(tool, False):
                results[key]['count'] += 1
                results[key][tool] = row['corrected_p']

            else:
                results[key] = {
                    'chr': row['#chr'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['alt'],
                    'count': 1,
                    tool: row['corrected_p']
                    }

def VarScan_process(tool, input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        freq = row.genotype('TUMOR')['FREQ']
        t_vaf = round(float(freq[:-1]) / 100, 4)

        if t_vaf < T_VAF_CUTOFF or float(row.INFO['SPV']) > P_VALUE_CUTOFF:
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
            t_vaf = round(float(row.INFO['TA']) / (float(row.INFO['TA']) + float(row.INFO['TR'])), 4)
        except ZeroDivisionError:
            t_vaf = 0

        if t_vaf < T_VAF_CUTOFF or float(row.INFO['PR']) < MUSEQ_CUTOFF:
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
            t_vaf = round(float(row.genotype('TUMOR')[str(row.ALT[0]) + 'U'][0])/float(row.genotype('TUMOR')['DP']), 4)
        except ZeroDivisionError:
            t_vaf = 0

        if t_vaf < T_VAF_CUTOFF or int(row.INFO['QSS']) < STRELKA_CUTOFF:
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

def VarScan_indel_process(input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        try:
            n_vaf = row.genotype('NORMAL')['AD']/row.genotype('NORMAL')['DP']
        except ZeroDivisionError:
            n_vaf = 0

        try:
            t_vaf = row.genotype('TUMOR')['AD']/row.genotype('TUMOR')['DP']
        except ZeroDivisionError:
            t_vaf = 0

        if t_vaf < T_VAF_CUTOFF or float(row.INFO['SPV']) > P_VALUE_CUTOFF:
            continue

        if results.get(key, False) and not results[key].get('VarScan', False):
            results[key]['VarScan'] = row.INFO['SPV']
            results[key]['N_coverage'] = row.genotype('NORMAL')['DP']
            results[key]['N_ref'] = row.genotype('NORMAL')['RD']
            results[key]['N_alt'] = row.genotype('NORMAL')['AD']
            results[key]['N_vaf'] = n_vaf
            results[key]['T_coverage'] = row.genotype('TUMOR')['DP']
            results[key]['T_ref'] = row.genotype('TUMOR')['RD']
            results[key]['T_alt'] = row.genotype('TUMOR')['AD']
            results[key]['T_vaf'] = t_vaf

        else:
            for alt in row.ALT:
                results[key] = {
                    'chr': row.CHROM,
                    'pos': str(row.POS),
                    'ref': str(row.REF),
                    'alt': str(alt),
                    'VarScan': str(row.INFO['SPV']),
                    'N_coverage': row.genotype('NORMAL')['DP'],
                    'N_ref': row.genotype('NORMAL')['RD'],
                    'N_alt': row.genotype('NORMAL')['AD'],
                    'N_vaf': n_vaf,
                    'T_coverage': row.genotype('TUMOR')['DP'],
                    'T_ref': row.genotype('TUMOR')['RD'],
                    'T_alt': row.genotype('TUMOR')['AD'],
                    'T_vaf': t_vaf
                    }

def Strelka_indel_process(input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        try:
            n_vaf = float(row.genotype('NORMAL')['TIR'][0])/float(row.genotype('NORMAL')['DP'])
        except ZeroDivisionError:
            n_vaf = 0
        try:
            t_vaf = float(row.genotype('TUMOR')['TIR'][0])/float(row.genotype('TUMOR')['DP'])
        except ZeroDivisionError:
            t_vaf = 0

        if t_vaf < T_VAF_CUTOFF or row.INFO['QSI'] < STRELKA_CUTOFF

        if results.get(key, False) and not results[key].get('Strelka', False):
            results[key]['Strelka'] = row.INFO['QSI']

        else:
            for alt in row.ALT:
                results[key] = {
                    'chr': row.CHROM,
                    'pos': str(row.POS),
                    'ref': str(row.REF),
                    'alt': str(alt),
                    'Strelka': str(row.INFO['QSI']),
                    'N_coverage': row.genotype('NORMAL')['DP'],
                    'N_ref': row.genotype('NORMAL')['TAR'][0],
                    'N_alt': row.genotype('NORMAL')['TIR'][0],
                    'N_vaf': n_vaf,
                    'T_coverage': row.genotype('TUMOR')['DP'],
                    'T_ref': row.genotype('TUMOR')['TAR'][0],
                    'T_alt': row.genotype('TUMOR')['TIR'][0],
                    'T_vaf': t_vaf
                    }



