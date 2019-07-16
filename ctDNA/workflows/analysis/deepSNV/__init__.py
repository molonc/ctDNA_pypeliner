import pypeliner
import pypeliner.managed as mgd

def run_deepSNV(config, normal_bam, tumour_bam, output_file):
    workflow = pypeliner.workflow.Workflow()

    workflow.commandline(
        name='r_deepSNV',
        args=(
            'Rscript',
            config["r_script_dir"] + 'deepSNV_analyze.R',
            '--tumour',
            tumour_bam,
            '--normal',
            normal_bam,
            '--bed',
            config["bed_file"],
            '--quality',
            25,
            '--out',
            output_file,
            )
        )

    return workflow
