import yaml
import argparse
import pypeliner
from workflows.analysis import deepSNV_workflow

if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('tumour', help='tumour bam file')
	argparser.add_argument('normal', help='normal bam file')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())

	pyp = pypeliner.app.Pypeline(config=args)
	workflow = pypeliner.workflow.Workflow()

	config = yaml.safe_load(open(args['config'], 'r'))
	tumour_bam = args['tumour']
	normal_bam = args['normal']
	tumour_sample = args['tumour'].split("/")[-1].split(".")[0]
	normal_sample = args['normal'].split("/")[-1].split(".")[0]
	workflow.subworkflow(
		name="analyze_with_deepSNV",
		func=deepSNV_workflow.run_deepSNV,
		args=(
			config,
			tumour_sample,
			tumour_bam,
			normal_bam,
			normal_sample,
			)
		)

	pyp.run(workflow)	