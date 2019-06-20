import csv
import vcf
from collections import OrderedDict

def create_result_dict(deepSNV_out, LoLoPicker_out, VarScan_out, museq_out, strelka_out):
    return {
    'deepSNV': {
        'file': deepSNV_out,
        'process_function': DeepSNV_process
        },
    'LoLoPicker': {
        'file': LoLoPicker_out,
        'process_function': LoLoPicker_process
        },
    'VarScan': {
        'file': VarScan_out,
        'process_function': vcf_process
        },
    'MutationSeq': {
        'file': museq_out,
        'process_function': vcf_process
        },
    'Strelka': {
        'file': strelka_out,
        'process_function': vcf_process
        }
    }

def merge_results(tumour_results, output_file):
    results = {}
    for tumour, result_file in tumour_results.iteritems():
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

        for key, result in sorted_results.iteritems():
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

        for key, result in sorted_results.iteritems():
            writer.writerow(result)

def DeepSNV_process(tool, input_file, results):
    with open(input_file, 'rb') as inputs:
        reader = csv.DictReader(inputs, delimiter='\t')
        for row in reader:
            key = row['chr'] + ':' + row['pos']
            if results.get(key, False) and not results[key].get(tool, False):
                results[key]['count'] += 1
                results[key][tool] = 'True'

            else:
                results[key] = {
                    'chr': row['chr'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['var'],
                    'count': 1,
                    tool: 'True'
                    }

def LoLoPicker_process(tool, input_file, results):
    with open(input_file, 'rb') as inputs:
        reader = csv.DictReader(inputs, delimiter='\t')
        for row in reader:
            key = row['#chr'] + ':' + row['pos']
            if results.get(key, False) and not results[key].get(tool, False):
                results[key]['count'] += 1
                results[key][tool] = 'True'

            else:
                results[key] = {
                    'chr': row['#chr'],
                    'pos': row['pos'],
                    'ref': row['ref'],
                    'alt': row['alt'],
                    'count': 1,
                    tool: 'True'
                    }

def vcf_process(tool, input_file, results):
    reader = vcf.Reader(filename=input_file)
    for row in reader:
        key = row.CHROM + ':' + str(row.POS)
        if results.get(key, False) and not results[key].get(tool, False):
            results[key]['count'] += 1
            results[key][tool] = 'True'

        else:
            for alt in row.ALT:
                results[key] = {
                    'chr': row.CHROM,
                    'pos': str(row.POS),
                    'ref': str(row.REF),
                    'alt': str(alt),
                    'count': 1,
                    tool: 'True'
                    }

def log_patient_analysis(input_files, output_file):
    with open(output_file, "w+") as output:
        output.write('tumour_sample\tresult_file\n')
        for tumour_id, result_file in input_files.iteritems():
            output.write(tumour_id + "\t" + result_file + "\n")

