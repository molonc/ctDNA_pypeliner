import yaml
import argparse
import pypeliner
import pypeliner.managed as mgd
from workflows import alignment

if __name__ == '__main__':
	argparser = argparse.ArgumentParser()
	pypeliner.app.add_arguments(argparser)

	argparser.add_argument('fastq_1', help='fastq file of read 1')
	argparser.add_argument('fastq_2', help='fastq file of read 2')
	argparser.add_argument('config', help='Configuration filename')

	args = vars(argparser.parse_args())

	pyp = pypeliner.app.Pypeline(config=args)
	workflow = pypeliner.workflow.Workflow()

	config = yaml.safe_load(open(args['config'], 'r'))
	fastq_1 = args['fastq_1']
	fastq_2 = args['fastq_2']
	sample_id = fastq_1.split("/")[-1].split("_")[0]

	workflow.subworkflow(
		name="align_sample",
		func=alignment.align_sample,
		args=(
			config,
			mgd.InputFile(fastq_1),
			mgd.InputFile(fastq_2),
			sample_id
			)
		)
	pyp.run(workflow)
