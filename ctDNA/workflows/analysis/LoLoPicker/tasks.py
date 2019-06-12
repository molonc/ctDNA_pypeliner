import os
from pypeliner.commandline import execute

def LoLoPicker_somatic(config, args, tumour_bam, normal_bam, somatic_file):
	execute(
		'LoLoPicker_somatic.py',
		'-t',
		tumour_bam,
		'-n',
		normal_bam,
		'-r',
		config['reference_genome'],
		'-b',
		config['bed_file'],
		'-o',
		args['results_dir'],
		)

def LoLoPicker_control(config, args, sample_list, somatic_file, control_file):
	execute(
		'LoLoPicker_control.py',
		'-l',
		sample_list,
		'-r',
		config['reference_genome'],
		'-o',
		args['results_dir'],
		)

def LoLoPicker_stats(args, control_file, stats_file, reject_file):
	execute(
			'LoLoPicker_stats.py',
			'-o',
			args['results_dir'],
			'--method',
			'FDR',
			)
	
def make_sample_list(normal_bams, sample_list_outfile):
	with open(sample_list_outfile, 'w+') as outfile:
		for sample, bam in normal_bams.items():
			outfile.write(bam + '\t')
			outfile.write(sample + '\n')