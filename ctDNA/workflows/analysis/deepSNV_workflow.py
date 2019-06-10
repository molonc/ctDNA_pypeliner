import pypeliner
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
