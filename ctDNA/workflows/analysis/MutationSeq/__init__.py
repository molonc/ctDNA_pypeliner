import pypeliner
import pypeliner.managed as mgd
import tasks

def run_MutationSeq(config, normal_bam, tumour_bam, output_file):
	workflow = pypeliner.workflow.Workflow()

	workflow.transform(
		name='generate_intervals',
		func=tasks.generate_intervals,
		ret=mgd.OutputChunks('interval'),
		args=(
			mgd.InputFile(config['bed_file']),
			)
		)

	workflow.transform(
		name='rum_museq_paired',
		axes=('interval',),
		func=tasks.run_museq,
		args=(
			config,
			mgd.InputFile(normal_bam),
			mgd.InputFile(tumour_bam),
			mgd.InputInstance('interval'),
			mgd.TempOutputFile('museq.vcf', 'interval'),
			mgd.TempOutputFile('museq.log', 'interval'),
			)
		)

	workflow.transform(
		name='merge_vcfs',
		func=tasks.merge_vcfs,
		args=(
			mgd.TempInputFile('museq.vcf', 'interval'),
			mgd.OutputFile(output_file),
			)
		)

	return workflow

