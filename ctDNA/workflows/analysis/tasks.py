import csv
from collections import OrderedDict

def merge_results(tumour_results, output_file):
    results = {}
    for result_file in tumour_results.itervalues():
        with open(result_file, 'rb') as result:
            reader = csv.DictReader(result, delimiter='\t')

            for row in reader:
                key = row['chr'] + ':' + row['pos']
                if results.get(key, False):
                    results[key]['count'] = int(results[key]['count']) + int(row['count'])

                else:
                    results[key] = row

    with open(output_file, 'w+') as output:
        field_names = [
            'chr',
            'pos',
            'ref',
            'alt',
            'count'
            ]

        writer = csv.DictWriter(
            output,
            fieldnames=field_names,
            restval=".",
            extrasaction='ignore',
            delimiter="\t",
            )
        writer.writeheader()

        sorted_results = OrderedDict(sorted(results.iteritems(), key=lambda x: (-int(x[1]['count']), x[1]['chr'], x[1]['pos'])))

        for result in sorted_results.itervalues():
            writer.writerow(result)

def log_patient_analysis(input_files, output_file):
    with open(output_file, "w+") as output:
        output.write('tumour_sample\tresult_file\n')
        for tumour_id, result_file in input_files.iteritems():
            output.write(tumour_id + "\t" + result_file + "\n")

