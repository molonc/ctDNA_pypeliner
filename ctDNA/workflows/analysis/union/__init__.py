import os
import csv
import tasks
from collections import OrderedDict

def create_result_dict(deepSNV_out, VarScan_out, museq_out, strelka_out, LoLoPicker_out):
    return {
    'deepSNV': {
        'file': deepSNV_out,
        'process_function': tasks.DeepSNV_process
        },
    'VarScan': {
        'file': VarScan_out,
        'process_function': tasks.VarScan_process
        },
    'MutationSeq': {
        'file': museq_out,
        'process_function': tasks.museq_process
        },
    'Strelka': {
        'file': strelka_out,
        'process_function': tasks.strelka_process
        },
    'LoLoPicker':{
            'file': LoLoPicker_out,
            'process_function': tasks.LoLoPicker_process
        }
    }

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
            'N_alf',
            'T_coverage',
            'T_A',
            'T_C',
            'T_G',
            'T_T',
            'T_N',
            'T_alf',
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
                tasks.bam_readcount(config, "N", normal_bam, result, os.path.join(union_space, 'normal_count.txt'))
                tasks.bam_readcount(config, "T", tumour_bam, result, os.path.join(union_space, 'tumour_count.txt'))
                if (result['T_coverage'] >= 1000 and
                    result['N_coverage'] >= 1000 and
                    result['T_alf'] > 0.004 and
                    result['T_alf'] > result['N_alf']):
                    writer.writerow(result)
