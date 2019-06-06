import os
import pypeliner
import pypeliner.managed
import tasks

def align_samples(
        config,
        fastq1_inputs,
        fastq2_inputs,
        bam_outputs,
        outdir):
    samples = bam_outputs.keys()

    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(
        obj=mgd.OutputChunks('sample_id'),
        value=samples
    )

    workflow.subworkflow(
        name='align_samples',
        func=align_sample,
        axes=('sample_id',),
        args=(
            config,
            mgd.InputFile('input.r1.fastq.gz', 'sample_id', fnames=fastq1_inputs),
            mgd.InputFile('input.r2.fastq.gz', 'sample_id', fnames=fastq2_inputs),
            mgd.OutputFile('output.bam', 'sample_id', fnames=bam_outputs),
            mgd.InputInstance("sample_id"),
            outdir
        ),
    )

    return workflow

def align_sample(config, fastq_1, fastq_2, out_file, outdir, ids):
    ref_genome = pypeliner.managed.InputFile(config['ref_genome']['file'])

    read_group_config = config.get('read_group', {})

    markdups_metrics = os.path.join(outdir, 'markdups_metrics.pdf')
    samtools_flagstat = os.path.join(outdir, 'samtools_flagstat.txt')

    out_bai = out_file + '.bai'

    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='align_bwa_mem',
        ctx={'mem': 8, 'ncpus': config['threads'], 'walltime': '08:00'},
        func=tasks.align_bwa_mem,
        args=(
            pypeliner.managed.TempInputFile('read_1', 'split'),
            pypeliner.managed.TempInputFile('read_2', 'split'),
            ref_genome,
            pypeliner.managed.TempOutputFile('aligned.bam', 'split'),
            config['threads'],
        ),
        kwargs={'sample_id': ids[0], 'lane_id': ids[1], 'read_group_info': config['read_group_info']}
    )

    workflow.transform(
        name='sort',
        ctx={'mem': 4, 'ncpus': 1, 'walltime': '08:00'},
        func=tasks.sort,
        args=(
            pypeliner.managed.TempInputFile('aligned.bam', 'split'),
            pypeliner.managed.TempOutputFile('sorted.bam', 'split'),
        ),
    )

    # workflow.transform(
    #     name='markdups',
    #     ctx={'mem': 8, 'ncpus': 1, 'walltime': '24:00'},
    #     func=tasks.markdups,
    #     args=(
    #         pypeliner.managed.TempInputFile('merged.bam'),
    #         pypeliner.managed.OutputFile(out_file),
    #         pypeliner.managed.OutputFile(markdups_metrics),
    #         pypeliner.managed.TempSpace("temp_markdups"),
    #     ),
    #     kwargs={'mem': '8G'},
    # )

    workflow.commandline(
        name='index',
        ctx={'mem': 4, 'ncpus': 1, 'walltime': '08:00'},
        args=(
            'samtools',
            'index',
            pypeliner.managed.InputFile(out_file),
            pypeliner.managed.OutputFile(out_bai)
        )
    )

    workflow.commandline(
        name='flagstat',
        ctx={'mem': 4, 'ncpus': 1, 'walltime': '08:00'},
        args=(
            'samtools',
            'flagstat',
            pypeliner.managed.InputFile(out_file),
            '>',
            pypeliner.managed.OutputFile(samtools_flagstat)
        )
    )

    return workflow