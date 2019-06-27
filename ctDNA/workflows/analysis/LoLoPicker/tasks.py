import os
from shutil import copyfile
from pypeliner.commandline import execute

def split_bed(bed_file):
    regions = []

    with open(bed_file, "rb") as bed:
        for region in bed.read().splitlines():
            regions.append(region.replace("\t", "-"))

    return regions

def create_axes_beds(region, bed_file):
    with open(bed_file, "wb") as bed:
        bed.write(region.replace("-", "\t"))
        bed.write("\n")

def make_sample_list(args, sample_list_outfile):
    normal_bams = args['normal_bams']

    with open(sample_list_outfile, 'w+') as outfile:
        for sample, bam in normal_bams.items():
            outfile.write(bam + '\t')
            outfile.write(sample + '\n')

def LoLoPicker_somatic(config, tumour_bam, normal_bam, region_bed, temp_dir, somatic_file):
    os.makedirs(temp_dir)

    execute(
        'LoLoPicker_somatic.py',
        '-t',
        tumour_bam,
        '-n',
        normal_bam,
        '-r',
        config['reference_genome'],
        '-b',
        region_bed,
        '-o',
        temp_dir,
        )

    copyfile(os.path.join(temp_dir, 'raw_somatic_varants.txt'), somatic_file)

def LoLoPicker_control(config, sample_list, temp_dir, somatic_file, control_file):
    os.makedirs(temp_dir)
    copyfile(somatic_file, os.path.join(temp_dir, 'raw_somatic_varants.txt'))

    execute(
        'LoLoPicker_control.py',
        '-l',
        sample_list,
        '-r',
        config['reference_genome'],
        '-o',
        temp_dir,
        )

    copyfile(os.path.join(temp_dir, 'control_stats.txt'), control_file)

def LoLoPicker_stats(temp_dir, somatic_file, control_file, output_file):
    os.makedirs(temp_dir)
    copyfile(somatic_file, os.path.join(temp_dir, 'raw_somatic_varants.txt'))
    copyfile(control_file, os.path.join(temp_dir, 'control_stats.txt'))

    execute(
            'LoLoPicker_stats.py',
            '-o',
            temp_dir,
            '--method',
            'FDR',
            )

    copyfile(os.path.join(temp_dir, 'stats_calls.txt'), output_file)

def _get_header(stats_calls):
    for line in stats_calls:
        if line.startswith("#"):
            return line

def merge_LoLoPicker(stats_calls_files, output_file):
    with open(output_file, "wb") as outfile:
        header = None

        for stats_calls_file in stats_calls_files.itervalues():
            with open(stats_calls_file, "rb") as stats_calls:
                if not header:
                    header = _get_header(stats_calls)

                    outfile.write(header)

                else:
                    _get_header(stats_calls)

                for line in stats_calls:
                    outfile.write(line)
