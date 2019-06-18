import os
from shutil import copyfile, copyfileobj
import gzip
from pypeliner.commandline import execute

def configure_bed(workspace, bed_file, output_bed, output_index):
    os.makedirs(workspace)

    execute(
        'bgzip',
        '-c',
        bed_file,
        '>',
        os.path.join(workspace, 'bed.gz'),
        )

    execute(
        'tabix',
        '-p',
        'bed',
        os.path.join(workspace, 'bed.gz'),
        )

    copyfile(os.path.join(workspace, 'bed.gz'), output_bed)
    copyfile(os.path.join(workspace, 'bed.gz.tbi'), output_index)

def run_strelka(config, normal_bam, tumour_bam, bed_file, bed_index, workspace, output_file):
    os.makedirs(workspace)

    copyfile(bed_file, os.path.join(workspace, 'bed.gz'))
    copyfile(bed_index, os.path.join(workspace, 'bed.gz.tbi'))

    execute(
        'configureStrelkaSomaticWorkflow.py',
        '--normalBam=' + normal_bam,
        '--tumorBam=' + tumour_bam,
        '--referenceFasta=' + config["reference_genome"],
        '--callRegions=' + os.path.join(workspace, 'bed.gz'),
        '--runDir=' + workspace,
        )

    execute(
        os.path.join(workspace, 'runWorkflow.py'),
        '--mode',
        'local',
        '--jobs',
        '8',
        '--quiet',
        )

    result_file = os.path.join(workspace, 'results/variants', 'somatic.snvs.vcf.gz')

    with gzip.open(result_file, 'rb') as results, open(output_file, 'wb') as output:
        copyfileobj(results, output)
