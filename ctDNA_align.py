import os
import shutil
import yaml
import glob
import argparse
import pypeliner
from pypeliner.commandline import execute
import pypeliner.workflow
import pypeliner.managed

def func_setup(fastq_dir):
    #import statements are within func during pypeliner dev. Once pypeliners is resolved correctly, remove import and place into global
    import os
    import shutil
    import glob

    alignment_list=glob.glob(fastq_dir+"*.fastq.gz")
    for align in alignment_list:
        sample_id=align.split("/")[-1].split("_")[0]
        print(sample_id)


if __name__ == '__main__':

    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('fastq_dir', help='dir of fastq' )
    # argparser.add_argument('config', help='Configuration Filename')

    args = vars(argparser.parse_args())
    
    #config = yaml.safe_load(open(args['config']).read())
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name='initialize',
        func=func_setup,
        args=(
            pypeliner.managed.InputFile(args["fastq_dir"]),
            #pypeliner.managed.OutputFile("fastq_files")
        )
	)

 #    workflow.commandline(
 #        name='align',
 #        args=(
 #            config['bwa'],
 #            'mem',
 #            config['reference'],
 #            config['fastq1'],
 #            config['fastq2'],
 #            '>',
 #            config['outfile'],
 #        )
	# )

    pyp.run(workflow)