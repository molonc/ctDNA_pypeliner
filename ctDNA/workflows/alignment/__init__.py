import pypeliner
import pypeliner.managed as mgd
import tasks

def align_sample(config, fastq_1, fastq_2, sample_id, out_bam, out_bai):
    workflow = pypeliner.workflow.Workflow()

    if config['umi_trim']:
        workflow.transform(
            name='trim_fastq',
            ctx={'mem': 8, 'ncpus': 1, 'walltime': '08:00'},
            func=tasks.trim_fastq,
            args=(
                mgd.InputFile(fastq_1),
                mgd.InputFile(fastq_2),
                mgd.TempSpace("trim_space"),
                mgd.TempOutputFile("fastq_1_trimmed.fastq"),
                mgd.TempOutputFile("fastq_2_trimmed.fastq"),
                )
            )
    else:
        workflow.transform(
            name='no_trim_fastq',
            ctx={'mem': 8, 'ncpus': 1, 'walltime': '08:00'},
            func=tasks.no_trim_fastq,
            args=(
                mgd.InputFile(fastq_1),
                mgd.InputFile(fastq_2),
                mgd.TempOutputFile("fastq_1_trimmed.fastq"),
                mgd.TempOutputFile("fastq_2_trimmed.fastq"),
                )
            )

    workflow.transform(
        name='fastq_to_sam',
        ctx={'mem': 8, 'ncpus': 1, 'walltime': '08:00'},
        func=tasks.fastq_to_sam,
        args=(
            mgd.InputFile(config["reference_genome"]),
            mgd.TempInputFile("fastq_1_trimmed.fastq"),
            mgd.TempInputFile("fastq_2_trimmed.fastq"),
            mgd.TempOutputFile('tmp.sam'),
            )
        )

    workflow.transform(
        name='sam_to_bam',
        func=tasks.sam_to_bam,
        args=(
            mgd.TempInputFile('tmp.sam'),
            mgd.TempOutputFile('tmp.bam'),
            )
        )

    workflow.transform(
        name='sort_bam',
        func=tasks.sort_bam,
        args=(
            mgd.TempInputFile('tmp.bam'),
            mgd.OutputFile(out_bam),
            )
        )

    workflow.transform(
        name='index_bam',
        func=tasks.index_bam,
        args=(
            mgd.InputFile(out_bam),
            mgd.OutputFile(out_bai),
            )
        )

    return workflow
