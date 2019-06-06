import yaml

def get_values_from_input(yamldata, key):
	values = {str(sample): yamldata[sample][key] for sample in yamldata}
	return values

def get_value_from_file(yamlfile, key):
	data = load_yaml(yamlfile)
	return data[key]

def get_fastq_files(input_yaml):
	data = load_yaml(input_yaml)

	items = {}
	for cell_id, cell_data in data.iteritems():
		items[cell_id] = {}
		for lane, laneinfo in cell_data["fastqs"].iteritems():
			items[cell_id][lane] = {}
			items[cell_id][lane]['fastq_1'] = format_file_yaml(laneinfo['fastq_1'])
			items[cell_id][lane]['fastq_2'] = format_file_yaml(laneinfo['fastq_2'])
	return items

def load_yaml(path):
	try:
		with open(path) as infile:
			data = yaml.safe_load(infile)

	except IOError:
		raise Exception(
			'Unable to open file: {0}'.format(path))
	return data
