import os
import pysam
from pypeliner.commandline import execute
from ctDNA.utils import vcfutils

def run_museq(config, normal_bam, tumour_bam, interval, output_file, log_file):
    reference = 'reference:' + config['reference_genome']
    model = 'model:' + config['museq_deep_model']

    tumour = 'tumour:' + tumour_bam
    normal = 'normal:' + normal_bam

    cmd = [config['museq_python'], config['museq_classify']]

    cmd.extend([reference, model, tumour, normal])
    cmd.extend([
        '--verbose',
        '--deep',
        '--purity', '70',
        '--coverage', config['coverage_threshold'],
        '--threshold', config['museq_threshold'],
        '--buffer_size', '2G',
        '--mapq_threshold', config['map_q'],
        '--indl_threshold', '0.05',
        '--normal_variant', int(config['N_vaf_threshold']) * 100,
        '--tumour_variant', '2',
        '--baseq_threshold', config['base_q'],
        '--config', config['museq_config'],
        '--manifest', config['bed_file']
        ])
    cmd.extend([
        '--out', output_file,
        '--log', log_file,
        '--interval', interval
        ])
    execute(*cmd)

def merge_vcfs(inputs, output_file, tempdir):
    os.makedirs(tempdir)
    mergedfile = os.path.join(tempdir, 'merged.vcf')
    vcfutils.concatenate_vcf(inputs, mergedfile)
    vcfutils.sort_vcf(mergedfile, output_file)
