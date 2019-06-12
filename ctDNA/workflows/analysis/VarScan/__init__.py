import pypeliner
import pypeliner.managed as mgd
import tasks

def run_VarScan(config, tumour_sample, tumour_bam, snp_output_file, indel_output_file):
	workflow = pypeliner.workflow.Workflow()

	workflow.transform(
		name='generate_mpileup',
		func=tasks.generate_mpileup,
		args=(
			config,
			mgd.InputFile(tumour_bam),
			mgd.TempOutputFile("tmp.mpileup"),
			)
		)

	workflow.transform(
		name='run_varscan_snp',
		func=tasks.run_varscan_snp,
		args=(
			mgd.TempInputFile("tmp.mpileup"),
			mgd.OutputFile(snp_output_file),
			)
		)

	workflow.transform(
		name='run_varscan_indel',
		func=tasks.run_varscan_indel,
		args=(
			mgd.TempInputFile("tmp.mpileup"),
			mgd.OutputFile(indel_output_file),
			)
		)

	return workflow

