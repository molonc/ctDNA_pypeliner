import os
import yaml
import argparse
import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

def deepSNV(args):
	pyp = pypeliner.app.Pypeline(modules=(), config=args)
	workflow = pypeliner.workflow.Workflow()

	config = yaml.safe_load(open(args['config'], 'r'))
	sample_id = args['tumour'].split("/")[-1].split("-")[0]

	workflow.commandline(
		name='run_deepSNV',
		args=(
			'Rscript',
			config["r_script_dir"] + 'deepSNV_analyze.R',
			'--tumour',
			args["tumour"],
			'--normal',
			args["normal"],
			'--bed',
			config["bed_file"],
			'--quality',
			10,
			'--out',
			mgd.OutputFile(config["results_dir"] + '{}_deepSNV_out.tsv'.format(sample_id)),

			)
		)

	pyp.run(workflow)

if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('tumour', help='tumour bam file')
	argparser.add_argument('normal', help='normal bam file')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())
	deepSNV(args)