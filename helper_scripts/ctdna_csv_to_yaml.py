import yaml
import csv

def main():
	yaml_dict = {}
	with open('ctDNA.csv', 'r') as csv_file:
		csv_dict = csv.DictReader(csv_file)
		for row in csv_dict:
			pbc_id = row["PBC ID"]
			if pbc_id not in yaml_dict:
				pbc_dict = {}
				pbc_dict["tumour"] = {}
				pbc_dict["normal"] = {}
				yaml_dict[pbc_id] = pbc_dict

			if row["Sample status"] == "Normal":
				sample_id = get_sample_sample_id(row)
				yaml_dict[pbc_id]["normal"][sample_id] = get_sample_dict(row)

			else:
				sample_id = get_sample_sample_id(row)
				yaml_dict[pbc_id]["tumour"][sample_id] = get_sample_dict(row)

	with open("ctDNA_mapping_no_fastq.yaml", "w+") as yaml_file:
		yaml.dump(yaml_dict, yaml_file, default_flow_style=False)

def get_sample_sample_id(row):
	aliquot_id = row["Aliquot ID"].replace(" ", "-")
	alias = row["Alias"].split("-")[0]
	return "-".join([aliquot_id, alias])

def get_sample_dict(row):
	sample_dict = {}
	sample_dict["run"] = row["Run Number"]
	sample_dict["type"] = row["Type"]
	sample_dict["sample_status"] = row["Sample status"]
	sample_dict["fastq1"] = ""
	sample_dict["fastq2"] = ""
	return sample_dict

if __name__ == '__main__':
	main()