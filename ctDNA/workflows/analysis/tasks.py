import os
from pypeliner.commandline import execute
from shutil import copyfile
from ctDNA.utils import vcfutils

def merge_normal(config, input_bams, output_file, output_bai):
    normal_list = list(bam for bam in input_bams.itervalues())
    cmd = ['samtools', 'merge', '-f', output_file]
    cmd.extend(normal_list)
    execute(*cmd)

    execute(
        'samtools',
        'index',
        output_file,
        output_bai,
        )

def annotate_outputs(config, temp_space, input_file, output_txt):
    os.makedirs(temp_space)

    execute(
        os.path.join(config['annovar'], 'convert2annovar.pl'),
        '-format',
        'vcf4',
        '-allsample',
        '-withfreq',
        input_file,
        '>',
        os.path.join(temp_space, 'anno_in')
        )

    execute(
        os.path.join(config['annovar'], 'table_annovar.pl'),
        os.path.join(temp_space, 'anno_in'),
        config['annovar_humandb'],
        '-buildver',
        'hg19',
        '-out',
        os.path.join(temp_space, 'ANNO'),
        '-remove',
        '-protocol',
        'refGene,cytoBand',
        '-operation',
        'g,r',
        '-nastring',
        '.',
        '-polish',
        )

    copyfile(os.path.join(temp_space, 'ANNO.hg19_multianno.txt'), output_txt)

def vcf_annotate_outputs(config, temp_space, input_file, output_vcf):
    os.makedirs(temp_space)

    execute(
        os.path.join(config['annovar'], 'table_annovar.pl'),
        input_file,
        config['annovar_humandb'],
        '-buildver',
        'hg19',
        '-out',
        os.path.join(temp_space, 'ANNO'),
        '-remove',
        '-protocol',
        'refGene,cytoBand',
        '-operation',
        'g,r',
        '-nastring',
        '.',
        '-polish',
        '-vcfinput'
        )

    copyfile(os.path.join(temp_space, 'ANNO.hg19_multianno.vcf'), output_vcf)

def log_patient_analysis(snv_tsv_files, indel_tsv_files, snv_txt_files, indel_txt_files, snv_vcf_files, indel_vcf_files, output_file):
    with open(output_file, "w+") as output:
        output.write('tumour_sample\tresult_file\n')
        for tumour_id, snv_tsv_file in snv_tsv_files.iteritems():
            output.write(tumour_id + "_snv_tsv\t" + snv_tsv_file + "\n")

        for tumour_id, indel_tsv_file in indel_tsv_files.iteritems():
            output.write(tumour_id + "_indel_tsv\t" + indel_tsv_file + "\n")

        for tumour_id, snv_txt_file in snv_txt_files.iteritems():
            output.write(tumour_id + "_snv_txt\t" + snv_txt_file + "\n")

        for tumour_id, indel_txt_file in indel_txt_files.iteritems():
            output.write(tumour_id + "_indel_txt\t" + indel_txt_file + "\n")

        for tumour_id, snv_vcf_file in snv_vcf_files.iteritems():
            output.write(tumour_id + "_snv_vcf\t" + snv_vcf_file + "\n")

        for tumour_id, indel_vcf_file in indel_vcf_files.iteritems():
            output.write(tumour_id + "_indel_vcf\t" + indel_vcf_file + "\n")
