import pypeliner
import pypeliner.managed as mgd
import deepSNV
import LoLoPicker
from collections import defaultdict

def partition_on_tumour(config, tumour_samples, normal_samples):
	workflow = pypeliner.workflow.Workflow()

	normal_bams = {str(sample): config["bam_directory"] + str(sample) + ".sorted.bam" for sample in normal_samples}
	tumour_bams = {str(sample): config["bam_directory"] + str(sample) + ".sorted.bam" for sample in tumour_samples}

	args = {
		'normal_samples': normal_samples,
		'normal_bams': normal_bams,
		'tumour_samples': tumour_samples,
		'tumour_bams': tumour_bams,
		'results_dir': config['results_dir']
		}


	workflow.setobj(obj=mgd.OutputChunks('tumour_id',), value=args['tumour_samples'])

	workflow.subworkflow(
		name='analyze_tumour',
		func=analyze_tumour,
		axes=('tumour_id',),
		args=(
			config,
			args,
			mgd.InputInstance('tumour_id'),
			mgd.InputFile('tumour_bam', 'tumour_id', fnames=args['tumour_bams']),
			)
		)

	return workflow

def analyze_tumour(config, args, tumour_sample, tumour_bam):
	workflow = pypeliner.workflow.Workflow()

	tumour_args = args

	tumour_args['results_dir'] = args['results_dir'] + tumour_sample + '/'

	workflow.commandline(
		name='make_tumour_results_directory',
		args=(
			'mkdir',
			'-p',
			tumour_args['results_dir']
			)
		)

	workflow.setobj(obj=mgd.OutputChunks('normal_id',), value=args['normal_samples'])

	workflow.subworkflow(
		name='analyze_tumour_normal',
		func=analyze_tumour_normal,
		axes=('normal_id',),
		args=(
			config,
			tumour_args,
			mgd.InputInstance('normal_id'),
			mgd.InputFile('normal_bam', 'normal_id', fnames=args['normal_bams']),
			tumour_sample,
			mgd.InputFile(tumour_bam),
			)
		)

	return workflow

def analyze_tumour_normal(config, args, normal_sample, normal_bam, tumour_sample, tumour_bam):
	workflow = pypeliner.workflow.Workflow()

	tumour_normal_args = args

	tumour_normal_args['results_dir'] = args['results_dir'] + '{}_{}/'.format(normal_sample, tumour_sample)

	workflow.commandline(
		name='make_tumour_normal_results_directory',
		args=(
			'mkdir',
			'-p',
			tumour_normal_args['results_dir']))

	workflow.subworkflow(
		name='run_deepSNV',
		func=deepSNV.run_deepSNV,
		args=(
			config,
			normal_sample,
			mgd.InputFile(normal_bam),
			tumour_sample,
			mgd.InputFile(tumour_bam),
			mgd.OutputFile(tumour_normal_args['results_dir'] + 'deepSNV_out.tsv')
			)
		)

	workflow.subworkflow(
		name='run_LoLoPicker',
		func=LoLoPicker.run_LoLoPicker,
		args=(
			config,
			tumour_normal_args,
			normal_sample,
			mgd.InputFile(normal_bam),
			tumour_sample,
			mgd.InputFile(tumour_bam),
			mgd.OutputFile(tumour_normal_args['results_dir'] + 'LoLoPicker_out.tsv'),
			)
		)

	return workflow
