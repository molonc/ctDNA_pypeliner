import os
import argparse
import pypeliner
import pypeliner.managed as mgd
from workflows import alignment
from workflows import analysis
from utils import helpers

def patient_workflow(config, patient_id, patient_input, output_file):
    workflow = pypeliner.workflow.Workflow()

    patient_bam_dir = config["bam_directory"] + patient_id
    patient_result_dir = config["results_dir"] + patient_id

    helpers.makedirs(patient_bam_dir)
    helpers.makedirs(patient_result_dir)

    input_args = helpers.create_input_args(patient_input, patient_bam_dir)

    workflow.setobj(obj=mgd.OutputChunks('sample_id',), value=input_args['all_samples'])

    workflow.subworkflow(
        name='align_samples',
        func=alignment.align_sample,
        axes=('sample_id',),
        args=(
            config,
            mgd.InputFile('fastq_1', 'sample_id', fnames=input_args['fastqs_r1']),
            mgd.InputFile('fastq_2', 'sample_id', fnames=input_args['fastqs_r2']),
            mgd.InputInstance('sample_id'),
            mgd.OutputFile('sample.bam', 'sample_id', fnames=input_args['all_bams']),
            mgd.OutputFile('sample.bam.bai', 'sample_id', fnames=input_args['all_bais']),
            )
        )

    workflow.subworkflow(
        name='run_analyses',
        func=analysis.partition_tumour,
        args=(
            config,
            input_args,
            patient_id,
            patient_result_dir,
            mgd.InputFile('sample.bam', 'sample_id', fnames=input_args['all_bams'], axes_origin=[]),
            mgd.InputFile('sample.bam.bai', 'sample_id', fnames=input_args['all_bais'], axes_origin=[]),
            mgd.OutputFile(output_file),
            )
        )

    return workflow

def ctDNA_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config'])
    for arg, value in args.iteritems():
        config[arg] = value

    helpers.makedirs(config["bam_directory"])

    helpers.makedirs(config["results_dir"])

    inputs = helpers.load_yaml(args['input_yaml'])
    patients = inputs.keys()

    workflow.setobj(obj=mgd.OutputChunks('patient_id',), value=patients)

    workflow.transform(
        name='get_input_by_patient',
        func=helpers.get_input_by_patient,
        ret=mgd.TempOutputObj('patient_input', 'patient_id'),
        axes=('patient_id',),
        args=(
            inputs,
            mgd.InputInstance('patient_id'),
            )
        )

    workflow.subworkflow(
        name='patient_workflow',
        func=patient_workflow,
        axes=('patient_id',),
        args=(
            config,
            mgd.InputInstance('patient_id'),
            mgd.TempInputObj('patient_input', 'patient_id'),
            mgd.OutputFile(os.path.join(config['results_dir'], '{patient_id}.log'), 'patient_id'),
            )
        )

    pyp.run(workflow)

def main():
    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument(
        '--input_yaml',
        required=True,
        help='input filename'
        )
    argparser.add_argument(
        '--config',
        required=True,
        help='Configuration filename'
        )
    argparser.add_argument(
        '--map_q',
        required=False,
        help='Minimum mapping quality',
        default=25
        )
    argparser.add_argument(
        '--base_q',
        required=False,
        help='Minimum base quality',
        default=15
        )
    argparser.add_argument(
        '--p_threshold',
        required=False,
        help='Maximum p_value (significance level)',
        default=0.0005
        )
    argparser.add_argument(
        '--museq_threshold',
        required=False,
        help='Minumum MutationSeq score',
        default=0.65
        )
    argparser.add_argument(
        '--N_vaf_threshold',
        required=False,
        help='Maximum normal variant allele frequency',
        default=0.6
        )
    argparser.add_argument(
        '--T_vaf_threshold',
        required=False,
        help='Minimum tumour variant allele frequency',
        default=0.004
        )
    argparser.add_argument(
        '--coverage_threshold',
        required=False,
        help='Minimum coverage',
        default=1000
        )
    argparser.add_argument(
        '--umi_trim',
        default=False,
        action='store_true'
        )

    args = vars(argparser.parse_args())
    ctDNA_workflow(args)


if __name__ == '__main__':
    main()