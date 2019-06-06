import os
import yaml
import argparse
import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd
import helpers
import workflows.alignment as alignment

def alignment_workflow(args):
	pyp = pypeliner.app.Pypeline(modules=(), config=args)
	workflow = pypeliner.workflow.Workflow()

	config = helpers.load_yaml(args['config'])
	inputs = helpers.load_yaml(args['input_yaml'])

	fastqs_r1 = helpers.get_values_from_input(inputs, 'fastq1')
	fastqs_r2 = helpers.get_values_from_input(inputs, 'fastq2')
	
	workflow.subworkflow(
		name="align_samples",
		func=alignment.align_samples,
		args=(
			config,
			fastqs_r2,
			fastqs_r1
			)
		)
	pyp.run(workflow)

if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('input_yaml', help='input filename')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())
	alignment_workflow(args)