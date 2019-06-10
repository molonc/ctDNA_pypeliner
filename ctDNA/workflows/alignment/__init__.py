import pypeliner
import pypeliner.managed as mgd
import tasks

def align_samples(config, fastq1_inputs, fastq2_inputs):
	samples = fastq1_inputs.keys()
	workflow = pypeliner.workflow.Workflow()

	workflow.setobj(obj=mgd.OutputChunks('sample_id',), value=samples)

	workflow.subworkflow(
		name='align_samples',
		func=align_sample,
		axes=('sample_id',),
		args=(
			config, 
			mgd.InputFile('fastq_1', 'sample_id', fnames=fastq1_inputs),
			mgd.InputFile('fastq_2', 'sample_id', fnames=fastq2_inputs),
			mgd.InputInstance('sample_id'),
			),
		)

	return workflow

def align_sample(config, fastq_1, fastq_2, sample_id):
	workflow = pypeliner.workflow.Workflow()

	workflow.transform(
		name='fastq_to_sam',
		func=tasks.fastq_to_sam,
		args=(
			config["reference_genome"],
			mgd.InputFile(fastq_1),
			mgd.InputFile(fastq_2),
			mgd.TempOutputFile('tmp.sam'),
			)
		)

	workflow.transform(
		name='sam_to_bam',
		func=tasks.sam_to_bam,
		args=(
			mgd.TempInputFile('tmp.sam'),
			mgd.TempOutputFile('tmp.bam'),
			)
		)

	workflow.transform(
		name='sort_bam',
		func=tasks.sort_bam,
		args=(
			mgd.TempInputFile('tmp.bam'),
			mgd.OutputFile(config["results_dir"] + '{}.sorted.bam'.format(sample_id))

			)
		)

	workflow.transform(
		name='index_bam',
		func=tasks.index_bam,
		args=(
			mgd.InputFile(config["results_dir"] + '{}.sorted.bam'.format(sample_id)),
			mgd.OutputFile(config["results_dir"] + '{}.sorted.bai'.format(sample_id)),
			)
		)

	return workflow
