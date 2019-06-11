from pypeliner.commandline import execute

def fastq_to_sam(ref_genome, fastq_1, fastq_2, output_sam):
	execute(
		'bwa',
		'mem',
		ref_genome,
		fastq_1,
		fastq_2,
		'>',
		output_sam
		)

def sam_to_bam(input_sam, output_bam):
	execute(
		'samtools',
		'view',
		'-bS',
		input_sam,
		'-o',
		output_bam
		)

def sort_bam(input_bam, output_sorted_bam):
	execute(
		'samtools',
		'sort',
		input_bam,
		'-o', 
		output_sorted_bam
		)

def index_bam(input_bam, output_bai):
	execute(
		'samtools', 
		'index',
		input_bam,
		output_bai
		)