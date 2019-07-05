import os
from shutil import copyfile
from pypeliner.commandline import execute
from ctDNA.utils import vcfutils

def create_axes_beds(bed_file, region, output_bed_file):
    with open(bed_file, "rb") as bed, open(output_bed_file, "wb") as output_bed:
        for line in bed:
            if line.startswith(region + "\t"):
                output_bed.write(line)


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
        '--mappingquality',
        10,
        '--basequality',
        20,
        '--tumoralteredreads',
        2,
        '--normalalteredreads',
        25,
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
        '--mappingquality',
        10,
        '--basequality',
        20,
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
            'Bonferroni',
            )

    copyfile(os.path.join(temp_dir, 'stats_calls.txt'), output_file)

def _get_header(stats_calls):
    for line in stats_calls:
        if line.startswith("#"):
            return line

def merge_LoLoPicker(temp_dir, stats_calls_files, output_file):
    os.makedirs(temp_dir)
    mergedfile = os.path.join(temp_dir, "merged.tsv")

    with open(mergedfile, "wb") as merged:
        header = None

        for stats_calls_file in stats_calls_files.itervalues():
            with open(stats_calls_file, "rb") as stats_calls:
                if not header:
                    header = _get_header(stats_calls)

                    merged.write(header)

                else:
                    _get_header(stats_calls)

                for line in stats_calls:
                    merged.write(line)

    vcfutils.sort_vcf(mergedfile, output_file)
