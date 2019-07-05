from __future__ import division
import os
import csv
import vcf
from collections import OrderedDict
from pypeliner.commandline import execute

def create_result_dict(deepSNV_out, VarScan_out, museq_out, strelka_out, LoLoPicker_out):
    return {
    'deepSNV': {
        'file': deepSNV_out,
        'process_function': DeepSNV_process
        },
    'VarScan': {
        'file': VarScan_out,
        'process_function': VarScan_process
        },
    'MutationSeq': {
        'file': museq_out,
        'process_function': museq_process
        },
    'Strelka': {
        'file': strelka_out,
        'process_function': strelka_process
        },
    'LoLoPicker':{
            'file': LoLoPicker_out,
            'process_function': LoLoPicker_process
        }
    }

def merge_results(tumour_results, output_file):
    results = {}
    for result_file in tumour_results.itervalues():
        with open(result_file, 'rb') as result:
            reader = csv.DictReader(result, delimiter='\t')

            for row in reader:
                key = row['chr'] + ':' + row['pos']
                if results.get(key, False):
                    results[key]['count'] = int(results[key]['count']) + int(row['count'])

                else:
                    results[key] = row

    with open(output_file, 'w+') as output:
        field_names = [
            'chr',
            'pos',
            'ref',
            'alt',
            'count'
            ]

        writer = csv.DictWriter(
            output,
            fieldnames=field_names,
            restval=".",
            extrasaction='ignore',
            delimiter="\t",
            )
        writer.writeheader()

        sorted_results = OrderedDict(sorted(results.iteritems(), key=lambda x: (-int(x[1]['count']), x[1]['chr'], x[1]['pos'])))

        for result in sorted_results.itervalues():
            writer.writerow(result)

def union_results(config, normal_bam, tumour_bam, tool_results, union_space, output_file):
    os.makedirs(union_space)
    results = {}
    for tool, process in tool_results.iteritems():
        process['process_function'](tool, process['file'], results)

    with open(output_file, 'w+') as output:
        field_names = [
            'chr',
            'pos',
            'ref',
            'alt',
            'count',
            'deepSNV',
            'LoLoPicker',
            'VarScan',
            'MutationSeq',
            'Strelka',
            'N_coverage',
            'N_A',
            'N_C',
            'N_G',
            'N_T',
            'N_N',
            'T_coverage',
            'T_A',
            'T_C',
            'T_G',
            'T_T',
            'T_N'
            ]

        writer = csv.DictWriter(
            output,
            fieldnames=field_names,
            restval=".",
            extrasaction='ignore',
            delimiter="\t",
            )
        writer.writeheader()

        sorted_results = OrderedDict(sorted(results.iteritems(), key=lambda x: (-x[1]['count'], x[1]['chr'], x[1]['pos'])))

        for result in sorted_results.itervalues():
            if result['count'] > 1:
                bam_readcount(config, "N", normal_bam, result, os.path.join(union_space, 'normal_count.txt'))
                bam_readcount(config, "T", tumour_bam, result, os.path.join(union_space, 'tumour_count.txt'))
                writer.writerow(result)

def bam_readcount(config, bam_type, bam, result, tmp_file):
    execute(
        'bam-readcount',
        '-d',
        20000,
        '-w',
        0,
        '--min-mapping-quality',
        20,
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
            result[bam_type + '_' + 'coverage'] = int(read_counts[3])
            for nucleotide_count in read_counts[5:10]:
                nucleotide_count_split = nucleotide_count.split(":")
                result[bam_type + '_' + nucleotide_count_split[0]] = int(nucleotide_count_split[1])

def DeepSNV_process(tool, input_file, results):
    with open(input_file, 'rb') as inputs:
        reader = csv.DictReader(inputs, delimiter='\t')
        for row in reader:
            key = row['chr'] + ':' + row['pos']
            try:
                taf = round(float(row['freq.var']), 4)
            except ZeroDivisionError:
                taf = 0

            if taf < 0.001:
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

        if taf < 0.001 or float(row.INFO['SPV']) > 0.05:
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

        if taf < 0.001:
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

        if taf < 0.001:
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

def log_patient_analysis(input_files, output_file):
    with open(output_file, "w+") as output:
        output.write('tumour_sample\tresult_file\n')
        for tumour_id, result_file in input_files.iteritems():
            output.write(tumour_id + "\t" + result_file + "\n")

