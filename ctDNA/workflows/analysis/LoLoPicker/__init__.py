import pypeliner
import pypeliner.managed as mgd
import tasks

def run_LoLoPicker(config, args, normal_sample, normal_bam, tumour_sample, tumour_bam, output_file):
	workflow = pypeliner.workflow.Workflow()

	workflow.commandline(
		name='LoLoPicker_somatic',
		func=tasks.LoLoPicker_somatic,
		args=(
			config, 
			args,
			mgd.InputFile(tumour_bam),
			mgd.InputFile(normal_bam),
			mgd.OutputFile(args['results_dir'] + "raw_somatic_variants.txt")
			)
		)

	workflow.transform(
		name='make_sample_list',
		func=tasks.make_sample_list,
		args=(
			args['normal_bams'],
			mgd.TempOutputFile('samplelist.txt'),
			)
		)

	workflow.commandline(
		name='LoLoPicker_control',
		func=tasks.LoLoPicker_control,
		args=(
			config,
			args,
			mgd.TempInputFile('samplelist.txt'),
			mgd.InputFile(args['results_dir'] + "raw_somatic_variants.txt"),
			mgd.OutputFile(args['results_dir'] + "control_stats.txt")
			)
		)

	workflow.commandline(
		name='LoLoPicker_stats',
		func=tasks.LoLoPicker_stats,
		args=(
			args['results_dir'],
			mgd.InputFile(args['results_dir'] + "control_stats.txt"),
			mgd.OutputFile(args['results_dir'] + "stats_calls.txt"),
			mgd.OutputFile(args['results_dir'] + "reject_calls.txt"),
			)
		)

	workflow.commandline(
		name='move_results',
		args=(
			'mv',
			mgd.InputFile(args['results_dir'] + 'stats_calls.txt'),
			output_file,
			)
		)

	return workflow
