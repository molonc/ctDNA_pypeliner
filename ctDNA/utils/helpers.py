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
