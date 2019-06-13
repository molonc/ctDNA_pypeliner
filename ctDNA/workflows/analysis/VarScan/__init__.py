import pypeliner
import pypeliner.managed as mgd
import tasks


def run_VarScan(config, normal_bam, tumour_bam, snp_output_file, indel_output_file):
	workflow = pypeliner.workflow.Workflow()

	workflow.transform(
		name='generate_normal_mpileup',
		func=tasks.generate_mpileup,
		args=(
			config,
			mgd.InputFile(normal_bam),
			mgd.TempOutputFile("normal.pileup"),
			)
		)

	workflow.transform(
		name='generate_tumour_mpileup',
		func=tasks.generate_mpileup,
		args=(
			config,
			mgd.InputFile(tumour_bam),
			mgd.TempOutputFile("tumour.pileup"),
			)
		)

	workflow.transform(
		name='run_varscan_somatic',
		func=tasks.run_varscan_somatic,
		args=(
			mgd.TempInputFile("normal.pileup"),
			mgd.TempInputFile("tumour.pileup"),
			mgd.OutputFile(snp_output_file),
			mgd.OutputFile(indel_output_file),
			)
		)

	return workflow

