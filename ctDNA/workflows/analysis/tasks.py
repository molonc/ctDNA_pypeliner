from pypeliner.commandline import execute

def merge_normal(config, input_bams, output_file, output_bai):
    normal_list = list(bam for bam in input_bams.itervalues())
    cmd = ['samtools', 'merge', output_file]
    cmd.extend(normal_list)
    execute(*cmd)

    execute(
        'samtools',
        'index',
        output_file,
        output_bai,
        )

def log_patient_analysis(input_files, output_file):
    with open(output_file, "w+") as output:
        output.write('tumour_sample\tresult_file\n')
        for tumour_id, result_file in input_files.iteritems():
            output.write(tumour_id + "\t" + result_file + "\n")
