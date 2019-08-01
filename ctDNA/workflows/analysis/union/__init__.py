import os
import csv
import vcf
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

def union_results(config, normal_bam, tumour_bam, tool_results, union_space, output_tsv, output_vcf):
    os.makedirs(union_space)
    results = {}
    for tool, process in tool_results.iteritems():
        process['process_function'](config, tool, process['file'], results)

    with open(output_tsv, 'wb') as tsv_output, open(output_vcf, 'wb') as vcf_output:
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
            'N_vaf',
            'T_coverage',
            'T_A',
            'T_C',
            'T_G',
            'T_T',
            'T_N',
            'T_vaf',
            ]

        tsv_writer = csv.DictWriter(
            tsv_output,
            fieldnames=field_names,
            restval='.',
            extrasaction='ignore',
            delimiter='\t',
            )

        tsv_writer.writeheader()


        vcf_template = vcf.Reader(filename=config['snv_vcf_template'])
        vcf_writer = vcf.Writer(vcf_output, vcf_template)

        sorted_results = OrderedDict(sorted(results.iteritems(), key=lambda x: (-x[1]['count'], x[1]['chr'], x[1]['pos'])))

        for result in sorted_results.itervalues():
            if result['count'] > 1:
                tasks.bam_readcount(config, "N", normal_bam, result, os.path.join(union_space, 'normal_count.txt'))
                tasks.bam_readcount(config, "T", tumour_bam, result, os.path.join(union_space, 'tumour_count.txt'))
                if (result['T_coverage'] >= config['coverage_threshold'] and
                    result['N_coverage'] >= config['coverage_threshold'] and
                    result['T_vaf'] > config['T_vaf_cutoff'] and
                    result['N_vaf'] < config['N_vaf_cutoff']):
                    tasks.write_snv_record(result, vcf_writer)
                    tsv_writer.writerow(result)

def union_indels(config, Strelka_in, VarScan_in,  output_tsv, output_vcf):
    results = {}

    tasks.Strelka_indel_process(config, Strelka_in, results)
    tasks.VarScan_indel_process(config, VarScan_in, results)

    with open(output_tsv, 'wb') as tsv_output, open(output_vcf, 'wb') as vcf_output:
        field_names = [
            'chr',
            'pos',
            'ref',
            'alt',
            'VarScan',
            'Strelka',
            'N_coverage',
            'N_ref',
            'N_alt',
            'N_vaf',
            'T_coverage',
            'T_ref',
            'T_alt',
            'T_vaf'
            ]

        tsv_writer = csv.DictWriter(
            tsv_output,
            fieldnames=field_names,
            restval='.',
            extrasaction='ignore',
            delimiter='\t',
            )

        tsv_writer.writeheader()

        vcf_template = vcf.Reader(filename=config['indel_vcf_template'])
        vcf_writer = vcf.Writer(vcf_output, vcf_template)

        sorted_results = OrderedDict(sorted(results.iteritems(), key=lambda x: (x[1]['chr'], x[1]['pos'])))

        for result in sorted_results.itervalues():
            tasks.write_indel_record(result, vcf_writer)
            tsv_writer.writerow(result)
