import os
import yaml
import argparse
import pypeliner
import pypeliner.workflow
import pypeliner.managed as mgd

def run_deepSNV(config, tumour_sample, tumour_bam, normal_sample, normal_bam):
	workflow = pypeliner.workflow.Workflow()

	workflow.commandline(
		name='r_deepSNV',
		args=(
			'Rscript',
			config["r_script_dir"] + 'deepSNV_analyze.R',
			'--tumour',
			mgd.InputFile(tumour_bam),
			'--normal',
			mgd.InputFile(normal_bam),
			'--bed',
			config["bed_file"],
			'--quality',
			10,
			'--out',
			mgd.OutputFile(config["results_dir"] + '{}-{}_deepSNV_out.tsv'.format(normal_sample, tumour_sample)),

			)
		)

	return workflow

if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('tumour', help='tumour bam file')
	argparser.add_argument('normal', help='normal bam file')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())

	pyp = pypeliner.app.Pypeline(modules=(), config=args)
	workflow = pypeliner.workflow.Workflow()

	config = yaml.safe_load(open(args['config'], 'r'))
	tumour_bam = args['tumour']
	normal_bam = args['normal']
	tumour_sample = args['tumour'].split("/")[-1].split(".")[0]
	normal_sample = args['normal'].split("/")[-1].split(".")[0]
	workflow.subworkflow(
		name="analyze_with_deepSNV",
		func=run_deepSNV,
		args=(
			config,
			tumour_sample,
			tumour_bam,
			normal_bam,
			normal_sample,
			)
		)

	pyp.run(workflow)	