from pypeliner.commandline import execute

def generate_mpileup(config, bam_file, output_file):
    execute(
        'samtools',
        'mpileup',
        '-B',
        '-Q',
        config['base_q'],
        '-C',
        '50',
        '-q',
        config['map_q'],
        '-d',
        '50000',
        '-f',
        config['reference_genome'],
        '-l',
        config['bed_file'],
        bam_file,
        '-o',
        output_file,
        )

def run_varscan_somatic(config, normal_mpileup, tumour_mpileup, snp_output_file, indel_outputfile):
    execute(
        'varscan',
        'somatic',
        normal_mpileup,
        tumour_mpileup,
        '--min-coverage',
        4,
        '--min-var-freq',
        config['T_vaf_threshold'],
        '--output-vcf',
        '1',
        '--output-snp',
        snp_output_file,
        '--output-indel',
        indel_outputfile,
        )
