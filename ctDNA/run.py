import argparse
import pypeliner
import pypeliner.managed as mgd
from workflows import alignment
from workflows import analysis
from utils import helpers

def create_input_args(input_file, config):
	import helpers

	inputs_raw = helpers.load_yaml(input_file)

	normal_samples = list(str(sample) for sample in inputs_raw["normal"])
	tumour_samples = list(str(sample) for sample in inputs_raw["tumour"])

	normal_bams = {str(sample): config["bam_directory"] + str(sample) + ".sorted.bam" for sample in normal_samples}
	tumour_bams = {str(sample): config["bam_directory"] + str(sample) + ".sorted.bam" for sample in tumour_samples}

	fastqs_r1 = helpers.get_fastq_files(inputs_raw, 'fastq1')
	fastqs_r2 = helpers.get_fastq_files(inputs_raw, 'fastq2')

	return {
		'fastqs_r1': fastqs_r1,
		'fastqs_r2': fastqs_r2,
		'normal_samples': normal_samples,
		'normal_bams': normal_bams,
		'tumour_samples': tumour_samples,
		'tumour_bams': tumour_bams,
		}


def ctDNA_workflow(args):
	pyp = pypeliner.app.Pypeline(config=args)
	workflow = pypeliner.workflow.Workflow()

	config = helpers.load_yaml(args['config'])

	workflow.transform(
		name="create_input_args",
		func=create_input_args,
		ret=mgd.TempOutputObj('inputs'),
		args=(
			mgd.InputFile(args['input_yaml']),
			config,
			)
		)

	workflow.subworkflow(
		name="align_samples",
		func=alignment.align_samples,
		args=(
			config,
			mgd.TempInputObj('inputs'),
			)
		)

	workflow.subworkflow(
		name="run_analyses",
		func=analysis.partition_on_tumour,
		args=(
			config,
			mgd.TempInputObj('inputs'),
			)
		)

	pyp.run(workflow)

def main():
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('input_yaml', help='input filename')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())
	ctDNA_workflow(args)


if __name__ == '__main__':
	main()