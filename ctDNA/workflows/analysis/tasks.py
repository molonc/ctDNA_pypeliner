import csv

def create_result_dict(deepSNV_out, LoLoPicker_out, VarScan_out, museq_out, strelka_out):
    return {
        'deepSNV': deepSNV_out,
        'LoLoPicker': LoLoPicker_out,
        'VarScan': VarScan_out,
        'MutationSeq': museq_out,
        'Strelka': strelka_out
    }

def union_results(tool_results, output_file):
    with open(output_file, 'w+') as output:
        for tool, result in tool_results.items():
            output.write(tool + result + '\n')

def mergestuff(tumour_results, output_file):
    with open(output_file, 'w+') as output:
        for tumour, result in tumour_results.items():
            output.write(tumour + result + '\n')