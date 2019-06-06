import os
import yaml
import argparse
import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

def align_samples(config, fastq1_inputs, fastq2_inputs):
	samples = fastq1_inputs.keys()
	workflow = pypeliner.workflow.Workflow()

	workflow.setobj(obj=mgd.OutputChunks('sample_id'), value=samples)

	workflow.subworkflow(
		name='align_samples',
		func=align_sample,
		axes=('sample_id',),
		args=(
			config, 
			mgd.InputFile('fastq_1', 'sample_id', fnames=fastq1_inputs),
			mgd.InputFile('fastq_2', 'sample_id', fnames=fastq2_inputs),
			mgd.InputInstance('sample_id')
			),
		)

	return workflow

def align_sample(config, fastq_1, fastq_2, sample_id):
	workflow = pypeliner.workflow.Workflow()

	workflow.commandline(
		name='fastq_to_sam',
		args=(
			'bwa',
			'mem',
			config["reference_genome"],
			mgd.InputFile(fastq_1),
			mgd.InputFile(fastq_2),
			'>',
			mgd.TempOutputFile(config["results_dir"] + '{}.sam'.format(sample_id)),

			)
		)

	workflow.commandline(
		name='sam_to_bam',
		args=(
			'samtools',
			'view',
			'-bS',
			mgd.TempInputFile(config["results_dir"] + '{}.sam'.format(sample_id)),
			'-o',
			mgd.TempOutputFile(config["results_dir"] + '{}.bam'.format(sample_id)),

			)
		)

	workflow.commandline(
		name='sort_bam',
		args=(
			'samtools',
			'sort',
			mgd.TempInputFile(config["results_dir"] + '{}.bam'.format(sample_id)),
			'-o',
			mgd.OutputFile(config["results_dir"] + '{}.sorted.bam'.format(sample_id))

			)
		)

	workflow.commandline(
		name='index_bam',
		args=(
			'samtools',
			'index',
			mgd.InputFile(config["results_dir"] + '{}.sorted.bam'.format(sample_id)),

			)
		)

	return workflow


if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('fastq_1', help='fastq file of read 1')
	argparser.add_argument('fastq_2', help='fastq file of read 2')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())

	pyp = pypeliner.app.Pypeline(modules=(), config=args)
	workflow = pypeliner.workflow.Workflow()

	config = yaml.safe_load(open(args['config'], 'r'))
	fastq_1 = args['fastq_1']
	fastq_2 = args['fastq_2']
	sample_id = fastq_1.split("/")[-1].split("_")[0]

	workflow.subworkflow(
		name="align_sample",
		func=align_sample,
		args=(
			config,
			mgd.InputFile(fastq_1),
			mgd.InputFile(fastq_2),
			sample_id
			)
		)
	pyp.run(workflow)

