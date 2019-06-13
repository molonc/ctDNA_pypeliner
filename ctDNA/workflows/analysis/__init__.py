import pypeliner
import pypeliner.managed as mgd
import deepSNV
import LoLoPicker
import VarScan
from collections import defaultdict

def partition_on_tumour(config, args):
	workflow = pypeliner.workflow.Workflow()

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

	tumour_results_dir = config['results_dir'] + tumour_sample + '/'

	workflow.commandline(
		name='make_tumour_results_directory',
		args=(
			'mkdir',
			'-p',
			tumour_results_dir
			)
		)

	workflow.setobj(obj=mgd.OutputChunks('normal_id',), value=args['normal_samples'])

	workflow.subworkflow(
		name='analyze_tumour_normal',
		func=analyze_tumour_normal,
		axes=('normal_id',),
		args=(
			config,
			args,
			tumour_results_dir,
			mgd.InputInstance('normal_id'),
			mgd.InputFile('normal_bam', 'normal_id', fnames=args['normal_bams']),
			tumour_sample,
			mgd.InputFile(tumour_bam),
			)
		)

	return workflow

def analyze_tumour_normal(config, args, results_dir, normal_sample, normal_bam, tumour_sample, tumour_bam):
	workflow = pypeliner.workflow.Workflow()

	matched_results_dir = results_dir + '{}_{}/'.format(normal_sample, tumour_sample)

	workflow.commandline(
		name='make_tumour_normal_results_directory',
		args=(
			'mkdir',
			'-p',
			matched_results_dir,
			)
		)

	workflow.subworkflow(
		name='run_deepSNV',
		func=deepSNV.run_deepSNV,
		args=(
			config,
			mgd.InputFile(normal_bam),
			mgd.InputFile(tumour_bam),
			mgd.OutputFile(matched_results_dir + 'deepSNV_out.tsv')
			)
		)

	workflow.subworkflow(
		name='run_LoLoPicker',
		func=LoLoPicker.run_LoLoPicker,
		args=(
			config,
			args,
			mgd.InputFile(normal_bam),
			mgd.InputFile(tumour_bam),
			mgd.OutputFile(matched_results_dir + 'LoLoPicker_out.tsv'),
			)
		)

	workflow.subworkflow(
		name='run_VarScan',
		func=VarScan.run_VarScan,
		args=(
			config,
			mgd.InputFile(normal_bam),
			mgd.InputFile(tumour_bam),
			mgd.OutputFile(matched_results_dir + 'VarScan_snp.vcf'),
			mgd.OutputFile(matched_results_dir + 'VarScan_indel.vcf'),
			))

	return workflow
