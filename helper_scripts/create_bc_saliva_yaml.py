import yaml

def load_yaml(path):
    try:
        with open(path) as infile:
            data = yaml.safe_load(infile)

    except IOError:
        raise Exception(
            'Unable to open file: {0}'.format(path))
    return data

def main():
    mapping = load_yaml("/home/pye/ctDNA_pypeliner/ctDNA_mapping.yaml")
    new_mapping = {}
    for patient, samples in mapping.iteritems():
        patient_mapping = {'normal': {}, 'tumour': {}}
        for sample, value in samples['normal'].iteritems():
            if value['type'] == 'buffy coat':
                patient_mapping['tumour'] = {sample: value}
            elif value['type'] == 'saliva':
                patient_mapping['normal'] = {sample:value}

        new_mapping[patient] = patient_mapping

    with open("/home/pye/ctDNA_pypeliner/test_files/bc_saliva.yaml", "w+") as yaml_file:
        yaml.dump(new_mapping, yaml_file, default_flow_style=False)



if __name__ == '__main__':
    main()