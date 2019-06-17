import argparse
import pypeliner
import pypeliner.managed as mgd
from workflows import alignment
from workflows import analysis
from utils import helpers

def patient_workflow(config, patient_input):
    workflow = pypeliner.workflow.Workflow()

    workflow.transform(
        name="create_input_args",
        func=helpers.create_input_args,
        ret=mgd.TempOutputObj('inputs'),
        args=(
            patient_input,
            config,
            )
        )

    workflow.subworkflow(
        name="align_samples",
        func=alignment.align_samples,
        args=(
            config,
            mgd.TempInputObj('inputs'),
            )
        )

    workflow.subworkflow(
        name="run_analyses",
        func=analysis.partition_on_tumour,
        args=(
            config,
            mgd.TempInputObj('inputs'),
            )
        )

    return workflow

def ctDNA_workflow(args):
    pyp = pypeliner.app.Pypeline(config=args)
    workflow = pypeliner.workflow.Workflow()

    config = helpers.load_yaml(args['config'])
    inputs = helpers.load_yaml(args['input_yaml'])
    patients = inputs.keys()

    workflow.setobj(obj=mgd.OutputChunks('patient_id',), value=patients)

    workflow.transform(
        name='get_input_by_patient',
        func=helpers.get_input_by_patient,
        ret=mgd.TempOutputObj('patient_input', 'patient_id'),
        axes=('patient_id',),
        args=(
            inputs,
            mgd.InputInstance('patient_id'),
            )
        )

    workflow.subworkflow(
        name='patient_workflow',
        func=patient_workflow,
        axes=('patient_id',),
        args=(
            config,
            mgd.TempInputObj('patient_input', 'patient_id')
            )
        )

    pyp.run(workflow)

def main():
    argparser = argparse.ArgumentParser()
    pypeliner.app.add_arguments(argparser)

    argparser.add_argument('input_yaml', help='input filename')
    argparser.add_argument('config', help='Configuration filename')

    args = vars(argparser.parse_args())
    ctDNA_workflow(args)


if __name__ == '__main__':
    main()