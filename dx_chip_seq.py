import argparse
import os
import subprocess

import dxpy


ENCODE_CHIP_SEQ_PROJECT = 'project-BJ7Kj200p2zX7qjB1j8001YK'

REPLICATES_FOLDER = '/replicates'
CONTROLS_FOLDER = '/controls'
PSEUDO_REPLICATES_FOLDER = '/pseudo_replicates'
SPP_FOLDER = '/spp'
IDR_FOLDER = '/idr'

def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Generate DNAnexus workflow for the ENCODE chip_seq pipeline.')

    ap.add_argument('-r', '--replicates',
                    help='Replicate files.',
                    nargs='+',
                    required=True)

    ap.add_argument('-c', '--controls',
                    help='Control files.',
                    nargs='+',
                    required=True)

    ap.add_argument('-p', '--project_name',
                    help='DNAnexus project name',
                    required=True)

    ap.add_argument('-d', '--duplicates_removed',
                    help='Flag indicating that duplicates have been removed',
                    action='store_true',
                    required=False)

    return ap.parse_args()

def find_reference_file_by_name(reference_name):
    '''Looks up a reference file by name in the project that holds common tools. From Joe Dale's code.'''

    found = dxpy.find_one_data_object(classname="file", name=reference_name,
                                      project=ENCODE_CHIP_SEQ_PROJECT,
                                      folder='/Reference Data',
                                      recurse=True,
                                      zero_ok=False, more_ok=False, return_handler=True)
    print "Resolved %s to %s" % (reference_name, found.get_id())
    return dxpy.dxlink(found)

def find_applet_by_name(applet_name):
    '''Looks up an applet by name in the project that holds tools.  From Joe Dale's code.'''

    found = dxpy.find_one_data_object(classname="applet", name=applet_name,
                                      project=ENCODE_CHIP_SEQ_PROJECT,
                                      zero_ok=False, more_ok=False, return_handler=True)
    print "Resolved %s to %s" % (applet_name, found.get_id())
    return found

def get_project(project_name):
    project = dxpy.find_projects(name=project_name, name_mode='glob', return_handler=True)

    project = [p for p in project]
    if len(project) < 1:
        project = dxpy.DXProject(dxpy.api.project_new({'name': project_name, 'summary': 'ChIP-Seq Pipeline'})['id'])
    elif len(project) > 1:
        print 'Found more than 1 project matching ' + project_name + '.'
        print 'Please provide a unique project!'
        sys.exit(1)
    else:
        project = project[0]

    return project

def populate_workflow(wf, replicates, controls, project_name, duplicates_removed):
    '''This function will populate the workflow for the ChIP-Seq Pipeline.'''
    if duplicates_removed:
        spp_app_name = 'spp_nodups'
    else:
        spp_app_name = 'spp'

    if len(controls) > 1:
        control_merge_input = {'input_bams': controls}
        stage_id = wf.add_stage(find_applet_by_name('merge_bams'), stage_input=control_merge_input, folder=CONTROLS_FOLDER)
        control_merge_output = dxpy.dxlink({'stage': stage_id, 'outputField': 'merged_bam'})
    else:
        control_merge_output = controls[0]

    replicate_merge_input = {'input_bams': replicates}
    stage_id = wf.add_stage(find_applet_by_name('merge_bams'), stage_input=replicate_merge_input, folder=REPLICATES_FOLDER)
    replicate_merge_output = dxpy.dxlink({'stage': stage_id, 'outputField': 'merged_bam'})

    pooled_replicate_v_control_spp_input = {'input_bam': replicate_merge_output, 'control_bam': control_merge_output}
    stage_id = wf.add_stage(find_applet_by_name(spp_app_name), stage_input=pooled_replicate_v_control_spp_input, folder=SPP_FOLDER)
    pooled_replicate_v_control_peaks = dxpy.dxlink({'stage': stage_id, 'outputField': 'peaks'})

    pooled_pseudo_replicate_input = {'input_bam': replicate_merge_output}
    stage_id = wf.add_stage(find_applet_by_name('pseudoreplicator'), stage_input=pooled_pseudo_replicate_input, folder=REPLICATES_FOLDER)
    pooled_pseudo_replicate_1 = dxpy.dxlink({'stage': stage_id, 'outputField': 'pseudoreplicate_bam1'})
    pooled_pseudo_replicate_2 = dxpy.dxlink({'stage': stage_id, 'outputField': 'pseudoreplicate_bam2'})

    pooled_pseudo_replicate_1_v_control_spp_input = {'input_bam': pooled_pseudo_replicate_1, 'control_bam': control_merge_output}
    stage_id = wf.add_stage(find_applet_by_name(spp_app_name), stage_input=pooled_pseudo_replicate_1_v_control_spp_input, folder=SPP_FOLDER)
    pooled_pseudo_replicate_peaks = [dxpy.dxlink({'stage': stage_id, 'outputField': 'peaks'})]

    pooled_pseudo_replicate_2_v_control_spp_input = {'input_bam': pooled_pseudo_replicate_2, 'control_bam': control_merge_output}
    stage_id = wf.add_stage(find_applet_by_name(spp_app_name), stage_input=pooled_pseudo_replicate_2_v_control_spp_input, folder=SPP_FOLDER)
    pooled_pseudo_replicate_peaks += [dxpy.dxlink({'stage': stage_id, 'outputField': 'peaks'})]

    replicates_v_controls_peaks = []
    pseudo_replicates_v_controls_peaks = []
    for replicate in replicates:
        pseudo_replicator_input = {'input_bam': replicate}
        stage_id = wf.add_stage(find_applet_by_name('pseudoreplicator'), stage_input=pseudo_replicator_input, folder=PSEUDO_REPLICATES_FOLDER)
        pseudo_replicate_1 = dxpy.dxlink({'stage': stage_id, 'outputField': 'pseudoreplicate_bam1'})
        pseudo_replicate_2 = dxpy.dxlink({'stage': stage_id, 'outputField': 'pseudoreplicate_bam2'})

        replicate_v_control_spp_input = {'input_bam': replicate, 'control_bam': control_merge_output}
        stage_id = wf.add_stage(find_applet_by_name(spp_app_name), stage_input=replicate_v_control_spp_input, folder=SPP_FOLDER)
        replicate_v_control_peaks = dxpy.dxlink({'stage': stage_id, 'outputField': 'peaks'})
        replicates_v_controls_peaks += [replicate_v_control_peaks]

        pseudo_replicate_1_v_control_spp_input = {'input_bam': pseudo_replicate_1, 'control_bam': control_merge_output}
        stage_id = wf.add_stage(find_applet_by_name(spp_app_name), stage_input=pseudo_replicate_1_v_control_spp_input, folder=SPP_FOLDER)
        pseudo_replicate_1_v_control_peaks = dxpy.dxlink({'stage': stage_id, 'outputField': 'peaks'})
        pseudo_replicates_v_controls_peaks += [pseudo_replicate_1_v_control_peaks]

        pseudo_replicate_2_v_control_spp_input = {'input_bam': pseudo_replicate_2, 'control_bam': control_merge_output}
        stage_id = wf.add_stage(find_applet_by_name(spp_app_name), stage_input=pseudo_replicate_2_v_control_spp_input, folder=SPP_FOLDER)
        pseudo_replicate_2_v_control_peaks = dxpy.dxlink({'stage': stage_id, 'outputField': 'peaks'})
        pseudo_replicates_v_controls_peaks += [pseudo_replicate_2_v_control_peaks]

    idr_input = {'replicate_peaks_files': replicates_v_controls_peaks,
                 'pseudo_replicate_peaks_files': pseudo_replicates_v_controls_peaks,
                 'pooled_replicate_peaks_file': pooled_replicate_v_control_peaks,
                 'pooled_pseudo_replicate_peaks_file': pooled_pseudo_replicate_peaks,
                 'replicate_peaks_threshold': 0.01,
                 'pseudo_replicate_peaks_threshold': 0.02,
                 'pooled_pseudo_replicate_peaks_threshold': 0.0025,
                 'output_prefix': project_name,
                 'ranking_measure': 'signal.value',
                 'genome_table_filename': 'genome_table.human.hg19.txt'}
    stage_id = wf.add_stage(find_applet_by_name('idr'), stage_input=idr_input, folder=IDR_FOLDER)


def copy_files(fids, project, folder):
    new_fids = []
    for fid in fids:
        (pid, fid) = fid.split(':')
        f = dxpy.DXFile(dxid=fid,project=pid)
        fn = f.describe()['name']
        found_file = dxpy.find_one_data_object(classname='file', project=project.get_id(), folder=folder, zero_ok=True, name=fn)

        if found_file is None:
            new_fids += [dxpy.dxlink(f.clone(project.get_id(), folder))]
        else:
            new_fids += [dxpy.dxlink(found_file)]

    return new_fids

if __name__ == '__main__':
    args = get_args()

    project = get_project(args.project_name)
    print 'Project: ' + project.describe()['name']
    project.new_folder(REPLICATES_FOLDER, True)
    project.new_folder(CONTROLS_FOLDER, True)
    replicates = copy_files(args.replicates, project, REPLICATES_FOLDER)
    controls = copy_files(args.controls, project, CONTROLS_FOLDER)

    # Now create a new workflow
    wf = dxpy.new_dxworkflow(title='dx_chip_seq',
                             name='ENCODE ChIP-Seq 2.0',
                             description='The ENCODE ChIP-Seq Pipeline 2.0',
                             project=project.get_id())
    populate_workflow(wf, replicates, controls, project.describe()['name'], args.duplicates_removed)

