from pypeliner.commandline import execute

def generate_mpileup(config, tumour_bam, output_file):
	execute(
		'samtools',
		'mpileup',
		'-B',
		'-Q',
		'20',
		'-C',
		'50',
		'-q',
		'20',
		'-d',
		'20000',
		'-f',
		config['reference_genome'],
		'-l',
		config['bed_file'],
		tumour_bam,
		'-o',
		output_file,
		)

def run_varscan_snp(input_mpileup, output_file):
	execute(
		'VarScan',
		'mpileup2snp',
		input_mpileup,
		'--min-coverage',
		'4',
		'--min-reads2',
		'2',
		'--min-avg-qual',
		'20',
		'--min-var-freq',
		'0.001',
		'--min-freq-for-hom',
		'90',
		'--output-vcf',
		'1',
		'>',
		output_file)

def run_varscan_indel(input_mpileup, output_file):
	execute(
		'VarScan',
		'mpileup2indel',
		input_mpileup,
		'--min-coverage',
		'4',
		'--min-reads2',
		'2',
		'--min-avg-qual',
		'20',
		'--min-var-freq',
		'0.001',
		'--min-freq-for-hom',
		'90',
		'--output-vcf',
		'1',
		'>',
		output_file)

