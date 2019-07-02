from __future__ import division
import csv
import vcf
from collections import OrderedDict

def create_result_dict(deepSNV_out, VarScan_out, museq_out, strelka_out, LoLoPicker_out=None):
    result_dict = {
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
        }
    }

    if LoLoPicker_out:
        result_dict['LoLoPicker'] ={
            'file': LoLoPicker_out,
            'process_function': LoLoPicker_process
        }

    return result_dict

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

def union_results(tool_results, output_file):
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
            'Strelka'
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
            writer.writerow(result)

def DeepSNV_process(tool, input_file, results):
    with open(input_file, 'rb') as inputs:
        reader = csv.DictReader(inputs, delimiter='\t')
        for row in reader:
            key = row['chr'] + ':' + row['pos']
            if results.get(key, False) and not results[key].get(tool, False):
                results[key]['count'] += 1
                results[key][tool] = row['freq.var']

            else:
                results[key] = {
                    'chr': row['chr'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['var'],
                    'count': 1,
                    tool: row['freq.var']
                    }

def LoLoPicker_process(tool, input_file, results):
    with open(input_file, 'rb') as inputs:
        reader = csv.DictReader(inputs, delimiter='\t')
        for row in reader:
            key = row['#chr'] + ':' + row['pos']
            if results.get(key, False) and not results[key].get(tool, False):
                results[key]['count'] += 1
                results[key][tool] = row['tumor_alf']

            else:
                results[key] = {
                    'chr': row['#chr'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['alt'],
                    'count': 1,
                    tool: row['tumor_alf']
                    }

def VarScan_process(tool, input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        freq = row.genotype('TUMOR')['FREQ']
        taf = float(freq[:-1]) / 100
        if results.get(key, False) and not results[key].get(tool, False):
            results[key]['count'] += 1
            results[key][tool] = taf

        else:
            for alt in row.ALT:
                results[key] = {
                    'chr': row.CHROM,
                    'pos': str(row.POS),
                    'ref': str(row.REF),
                    'alt': str(alt),
                    'count': 1,
                    tool: taf
                    }

def museq_process(tool, input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        taf = float(row.INFO['TA']) / (float(row.INFO['TA']) + float(row.INFO['TR']))
        if results.get(key, False) and not results[key].get(tool, False):
            results[key]['count'] += 1
            results[key][tool] = taf

        else:
            for alt in row.ALT:
                results[key] = {
                    'chr': row.CHROM,
                    'pos': str(row.POS),
                    'ref': str(row.REF),
                    'alt': str(alt),
                    'count': 1,
                    tool: taf
                    }

def strelka_process(tool, input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        taf = float(row.genotype('TUMOR')[str(row.ALT[0]) + 'U'][0])/float(row.genotype('TUMOR')['DP'])
        if results.get(key, False) and not results[key].get(tool, False):
            results[key]['count'] += 1
            results[key][tool] = taf

        else:
            for alt in row.ALT:
                results[key] = {
                    'chr': row.CHROM,
                    'pos': str(row.POS),
                    'ref': str(row.REF),
                    'alt': str(alt),
                    'count': 1,
                    tool: taf
                    }

def log_patient_analysis(input_files, output_file):
    with open(output_file, "w+") as output:
        output.write('tumour_sample\tresult_file\n')
        for tumour_id, result_file in input_files.iteritems():
            output.write(tumour_id + "\t" + result_file + "\n")

