import os
import yaml
import argparse
import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

def loloPicker(args):
	pyp = pypeliner.app.Pypeline(modules=(), config=args)
	workflow = pypeliner.workflow.Workflow()

	config = yaml.safe_load(open(args['config'], 'r'))
	sample_id = args['tumour'].split("/")[-1].split("-")[0]

	workflow.commandline(
		
		)

	pyp.run(workflow)

if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('tumour', help='tumour bam file')
	argparser.add_argument('normal', help='normal bam file')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())
	loloPicker(args)