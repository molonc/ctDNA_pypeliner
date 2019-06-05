import os
import yaml
import argparse
import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

def align():
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('fastq_one', help='fastq file of read 1')
	argparser.add_argument('fastq_two', help='fastq file of read 2')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())
	config = yaml.safe_load(open(args['config'], 'r'))
	sample_id = args['fastq_one'].split("/")[-1].split("_")[0]

	pyp = pypeliner.app.Pypeline(modules=(), config=args)
	workflow = pypeliner.workflow.Workflow()

	workflow.commandline(
		name='fastq_to_sam',
		args=(
			'bwa',
			'mem',
			config["reference_genome"],
			args["fastq_one"],
			args["fastq_two"],
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
			mgd.OutputFile(config["results_dir"] + '{}.sorted.bam'.format(sample_id)),

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

	pyp.run(workflow)


if __name__ == '__main__':
	align()