import pypeliner
import pypeliner.managed as mgd
import deepSNV_workflow as deepSNV

def run_multi(config, tumour_samples, normal_samples):
	workflow = pypeliner.workflow.Workflow()

	normal_bams = {str(sample): config["results_dir"] + str(sample) + ".sorted.bam" for sample in normal_samples}
	tumour_bams = {str(sample): config["results_dir"] + str(sample) + ".sorted.bam" for sample in tumour_samples}

	workflow.setobj(obj=mgd.OutputChunks('normal_id',), value=normal_samples)

	workflow.subworkflow(
		name='run_deepSNV_on_normal',
		func=partition_on_normal,
		axes=('normal_id',),
		args=(
			config,
			tumour_samples,
			tumour_bams,
			mgd.InputInstance('normal_id'),
			mgd.InputFile('normal_bam', 'normal_id', fnames=normal_bams),
			)
		)

	return workflow

def partition_on_normal(config, tumour_samples, tumour_bams, normal_sample, normal_bam):
	workflow = pypeliner.workflow.Workflow()

	workflow.setobj(obj=mgd.OutputChunks('tumour_id',), value=tumour_samples)

	workflow.subworkflow(
		name='run_deepSNV',
		func=deepSNV.run_deepSNV,
		axes=('tumour_id',),
		args=(
			config,
			mgd.InputInstance('tumour_id'),
			mgd.InputFile('tumour_bam', 'tumour_id', fnames=tumour_bams),
			normal_sample,
			mgd.InputFile(normal_bam),
			)
		)

	return workflow