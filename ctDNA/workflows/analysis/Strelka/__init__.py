import pypeliner
import pypeliner.managed as mgd
import tasks

def run_Strelka(config, normal_bam, tumour_bam, output_file):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='configure_bed',
        func=tasks.configure_bed,
        args=(
            mgd.TempSpace('bed_space'),
            mgd.InputFile(config['bed_file']),
            mgd.TempOutputFile('bed.gz'),
            mgd.TempOutputFile('bed.gz.tbi')
            )
        )

    workflow.transform(
        name='run_strelka',
        func=tasks.run_strelka,
        args=(
            config,
            mgd.InputFile(normal_bam),
            mgd.InputFile(tumour_bam),
            mgd.TempInputFile('bed.gz'),
            mgd.TempInputFile('bed.gz.tbi'),
            mgd.TempSpace('strelka_workspace'),
            mgd.OutputFile(output_file)
            )
        )

    return workflow
