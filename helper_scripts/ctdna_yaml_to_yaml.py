import yaml
import re
from os import listdir
from os.path import join

DIR_PATH = "/shahlab/archive/miyuen_tmp/ctDNA_data"
# DIR_PATH = "/Users/pye/projects/Data/run23/TNBC1194/"
def main():
    fastqs = [fastq for fastq in listdir(DIR_PATH) if fastq.endswith("fastq.gz")]
    yaml_dict = {}
    with open('ctDNA_mapping_no_fastq.yaml', 'r') as no_fastq_yaml_file:
    # with open('test.yaml', 'r') as no_fastq_yaml_file:
        yaml_dict = yaml.safe_load(no_fastq_yaml_file)
        for fastq in fastqs:
            for pbc_dict in yaml_dict.itervalues():
                for sample_type_dict in pbc_dict.itervalues():
                    for sample_dict in sample_type_dict.itervalues():
                        if fastq.startswith(sample_id):
                            if re.search("_R1_", fastq):
                                sample_dict["fastq1"] = join(DIR_PATH, fastq)
                            elif re.search("_R2_", fastq):
                                sample_dict["fastq2"] = join(DIR_PATH, fastq)

    with open("ctDNA_mapping.yaml", "w+") as yaml_file:
        yaml.dump(yaml_dict, yaml_file, default_flow_style=False)


if __name__ == '__main__':
    main()