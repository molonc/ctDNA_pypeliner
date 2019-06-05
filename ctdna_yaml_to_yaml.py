import yaml
import re
from os import listdir
from os.path import isfile, join

DIR_PATH = "/Users/pye/beast_share/lustre/archive/MiSeq/MiSeq_Analysis_Files/190514_M02348_0006_000000000-C8WC4/Data/Intensities/BaseCalls/"
# DIR_PATH = "/Users/pye/projects/Data/run23/TNBC1194/"
def main():
	fastqs = [fastq for fastq in listdir(DIR_PATH) if fastq.endswith("fastq.gz")]
	yaml_dict = {}
	# with open('ctDNA_mapping_no_fastq.yaml', 'r') as no_fastq_yaml_file:
	with open('test.yaml', 'r') as no_fastq_yaml_file:
		yaml_dict = yaml.safe_load(no_fastq_yaml_file)
		for fastq in fastqs:
			for pbc, pbc_dict in yaml_dict.iteritems():
				for sample_type, sample_type_dict in pbc_dict.iteritems():
					for sample_id, sample_dict in sample_type_dict.iteritems():
						if fastq.startswith(sample_id):
							if re.search("_R1_", fastq):
								sample_dict["fastq1"] = join(DIR_PATH, fastq)
							elif re.search("_R2_", fastq):
								sample_dict["fastq2"] = join(DIR_PATH, fastq)

	with open("ctDNA_mapping.yaml", "w+") as yaml_file:
		yaml.dump(yaml_dict, yaml_file, default_flow_style=False)


if __name__ == '__main__':
	main()