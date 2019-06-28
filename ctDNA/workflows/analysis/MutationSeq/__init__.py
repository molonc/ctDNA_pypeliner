import pypeliner
import pypeliner.managed as mgd
import tasks

def run_MutationSeq(config, normal_bam, tumour_bam, output_file):
    workflow = pypeliner.workflow.Workflow()

    workflow.setobj(obj=mgd.OutputChunks('interval',), value=list(map(str, range(1, 23) + ['X'])))

    workflow.transform(
        name='run_museq_paired',
        ctx={'mem': 8, 'ncpus': 1, 'walltime': '24:00'},
        axes=('interval',),
        func=tasks.run_museq,
        args=(
            config,
            mgd.InputFile(normal_bam),
            mgd.InputFile(tumour_bam),
            mgd.InputInstance('interval'),
            mgd.TempOutputFile('museq.vcf', 'interval'),
            mgd.TempOutputFile('museq.log', 'interval'),
            )
        )

    workflow.transform(
        name='merge_vcfs',
        func=tasks.merge_vcfs,
        args=(
            mgd.TempInputFile('museq.vcf', 'interval', axes_origin=[]),
            mgd.OutputFile(output_file),
            mgd.TempSpace('merge_vcf'),
            )
        )

    return workflow

