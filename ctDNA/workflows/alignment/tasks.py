import os
from shutil import copyfile
from pypeliner.commandline import execute

def no_trim_fastq(fastq_1, fastq_2, output_fastq_1, output_fastq_2):
    copyfile(fastq_1, output_fastq_1)
    copyfile(fastq_2, output_fastq_2)

def trim_fastq(fastq_1, fastq_2, temp_dir, output_fastq_1, output_fastq_2):
    os.makedirs(temp_dir)

    execute(
        'trim_galore',
        '--paired',
        '--output_dir',
        temp_dir,
        '--clip_R1',
        13,
        '--three_prime_clip_R1',
        13,
        '--clip_R2',
        13,
        '--three_prime_clip_R2',
        13,
        fastq_1,
        fastq_2
        )

    trimmed_fastq_1 = os.path.basename(fastq_1).replace(".fastq.gz", "_val_1.fq.gz")
    trimmed_fastq_2 = os.path.basename(fastq_2).replace(".fastq.gz", "_val_2.fq.gz")

    copyfile(os.path.join(temp_dir, trimmed_fastq_1), output_fastq_1)
    copyfile(os.path.join(temp_dir, trimmed_fastq_2), output_fastq_2)

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