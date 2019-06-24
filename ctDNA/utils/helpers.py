import os
import errno
import yaml

def get_values_from_input(yamldata, key):
    values = {str(sample): yamldata[sample][key] for sample in yamldata}
    return values

def get_value_from_file(yamlfile, key):
    data = load_yaml(yamlfile)
    return data[key]

def get_fastq_files(data, key):

    items = {}
    for cell_type, cell_data in data.iteritems():
        for sample, sample_info in cell_data.iteritems():
            items[sample] = sample_info[key]
    return items

def load_yaml(path):
    try:
        with open(path) as infile:
            data = yaml.safe_load(infile)

    except IOError:
        raise Exception(
            'Unable to open file: {0}'.format(path))
    return data

def create_input_args(patient_input, patient_bam_dir, no_lolo):
    normal_samples = list(str(sample) for sample in patient_input["normal"])
    tumour_samples = list(str(sample) for sample in patient_input["tumour"])

    normal_bams = {str(sample): patient_bam_dir + str(sample) + ".sorted.bam" for sample in normal_samples}
    tumour_bams = {str(sample): patient_bam_dir + str(sample) + ".sorted.bam" for sample in tumour_samples}

    all_samples = normal_samples + tumour_samples
    all_bams = dict(normal_bams.items() + tumour_bams.items())

    fastqs_r1 = get_fastq_files(patient_input, 'fastq1')
    fastqs_r2 = get_fastq_files(patient_input, 'fastq2')

    return {
        'fastqs_r1': fastqs_r1,
        'fastqs_r2': fastqs_r2,
        'normal_samples': normal_samples,
        'normal_bams': normal_bams,
        'tumour_samples': tumour_samples,
        'tumour_bams': tumour_bams,
        'all_samples': all_samples,
        'all_bams': all_bams,
        'run_LoLoPicker': not no_lolo
        }

def get_input_by_patient(inputs, patient_id):
    return inputs[patient_id]

def makedirs(dir_name):
    if not os.path.exists(dir_name):
        try:
            os.makedirs(dir_name)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
