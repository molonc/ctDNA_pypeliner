import pypeliner
import pypeliner.managed as mgd
import tasks

def run_LoLoPicker(config, args, normal_bam, tumour_bam, output_file):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='LoLoPicker_somatic',
        func=tasks.LoLoPicker_somatic,
        args=(
            config,
            mgd.InputFile(tumour_bam),
            mgd.InputFile(normal_bam),
            mgd.TempSpace('LoLoPicker_somatic_temp'),
            mgd.TempOutputFile("raw_somatic_variants.txt")
            )
        )

    workflow.transform(
        name='make_sample_list',
        func=tasks.make_sample_list,
        args=(
            args,
            mgd.TempOutputFile('samplelist.txt'),
            )
        )

    workflow.transform(
        name='LoLoPicker_control',
        func=tasks.LoLoPicker_control,
        args=(
            config,
            mgd.TempInputFile('samplelist.txt'),
            mgd.TempSpace('LoLoPicker_control_temp'),
            mgd.TempInputFile("raw_somatic_variants.txt"),
            mgd.TempOutputFile("control_stats.txt")
            )
        )

    workflow.transform(
        name='LoLoPicker_stats',
        func=tasks.LoLoPicker_stats,
        args=(
            mgd.TempSpace('LoLoPicker_stats_temp'),
            mgd.TempInputFile("raw_somatic_variants.txt"),
            mgd.TempInputFile("control_stats.txt"),
            mgd.OutputFile(output_file),
            )
        )

    return workflow
