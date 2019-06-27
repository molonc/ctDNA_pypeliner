import pypeliner
import pypeliner.managed as mgd
import tasks

def run_LoLoPicker(config, args, normal_bam, tumour_bam, output_file):
    workflow = pypeliner.workflow.Workflow()

    regions = tasks.split_bed(config["bed_file"])

    workflow.setobj(obj=mgd.OutputChunks('region',), value=regions)

    workflow.transform(
        name='create_axes_beds',
        axes=('region',),
        func=tasks.create_axes_beds,
        args=(
            mgd.InputInstance('region'),
            mgd.TempOutputFile('region.bed', 'region')
            )
        )

    workflow.transform(
        name='LoLoPicker_somatic',
        axes=('region',),
        func=tasks.LoLoPicker_somatic,
        args=(
            config,
            mgd.InputFile(tumour_bam),
            mgd.InputFile(normal_bam),
            mgd.TempInputFile('region.bed', 'region'),
            mgd.TempSpace('LoLoPicker_somatic_temp', 'region'),
            mgd.TempOutputFile("raw_somatic_varants.txt", 'region')
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
        axes=('region',),
        func=tasks.LoLoPicker_control,
        args=(
            config,
            mgd.TempInputFile('samplelist.txt'),
            mgd.TempSpace('LoLoPicker_control_temp', 'region'),
            mgd.TempInputFile("raw_somatic_varants.txt", 'region'),
            mgd.TempOutputFile("control_stats.txt", 'region')
            )
        )

    workflow.transform(
        name='LoLoPicker_stats',
        axes=('region',),
        func=tasks.LoLoPicker_stats,
        args=(
            mgd.TempSpace('LoLoPicker_stats_temp', 'region'),
            mgd.TempInputFile("raw_somatic_varants.txt", 'region'),
            mgd.TempInputFile("control_stats.txt", 'region'),
            mgd.TempOutputFile("stats_calls.txt", 'region'),
            )
        )

    workflow.transform(
        name='merge_LoLoPicker',
        func=tasks.merge_LoLoPicker,
        args=(
            mgd.TempInputFile("stats_calls.txt", 'region', axes_origin=[]),
            mgd.OutputFile(output_file)
            )
        )

    return workflow
