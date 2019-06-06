import os
import shutil
import yaml
import glob
import argparse
import pypeliner
from pypeliner.commandline import execute
import pypeliner.workflow
import pypeliner.managed

if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('fastq_dir', help='fastq dir')
    argparser.add_argument('out_dir', help='output dir')
    argparser.add_argument('config', help='Configuration Filename')

    args = vars(argparser.parse_args())
    
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = yaml.safe_load(open(args['config']).read())
    fastqs_r1 = helpers.get_values_from_input(args['fastq_dir'], 'fastq1')
    fastqs_r2 = helpers.get_values_from_input(args['fastq_dir'], 'fastq2')
    outputs = helpers.get_values_from_input(args['fastq_dir'], 'bam')
    outdir = args['out_dir']

    workflow.subworkflow(
        name="align_samples",
        func=alignment.align_samples,
        args=(
            config,
            fastqs_r1,
            fastqs_r2,
            outputs,
            outdir
        )
    )

    pyp.run(workflow)