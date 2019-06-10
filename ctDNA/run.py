import argparse
import pypeliner
import helpers
from workflows import alignment
from workflows import analysis

def ctDNA_workflow(args):
	pyp = pypeliner.app.Pypeline(config=args)
	workflow = pypeliner.workflow.Workflow()

	config = helpers.load_yaml(args['config'])
	inputs = helpers.load_yaml(args['input_yaml'])

	normal_samples = list(str(sample) for sample in inputs["normal"])
	tumour_samples = list(str(sample) for sample in inputs["tumour"])

	fastqs_r1 = helpers.get_fastq_files(inputs, 'fastq1')
	fastqs_r2 = helpers.get_fastq_files(inputs, 'fastq2')

	workflow.subworkflow(
		name="align_samples",
		func=alignment.align_samples,
		args=(
			config,
			fastqs_r1,
			fastqs_r2,
			)
		)

	workflow.subworkflow(
		name="run_anlyses",
		func=analysis.run_multi,
		args=(
			config, 
			tumour_samples,
			normal_samples,
			)
		)

	pyp.run(workflow)

if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('input_yaml', help='input filename')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())
	ctDNA_workflow(args)