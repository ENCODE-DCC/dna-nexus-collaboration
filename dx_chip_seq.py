import argparse
import os
import subprocess

import dxpy

REFERENCE_DATA = {'reference_genome':           'human_g1k_v37_decoy.fasta.gz',
                  'reference_genome_bwa_index': 'human_g1k_v37_decoy.fasta.bwa.indexed.tar.gz',
                  'reference_genome_fai':       'human_g1k_v37_decoy.fasta.fai',
                  'reference_genome_dict':      'human_g1k_v37_decoy.dict'}

ENCODE_CHIP_SEQ_PROJECT = ''

INPUT_FASTQS_FOLDER = '/input_fastqs'

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
                                      folder='/Applets',
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

def get_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Generate DNAnexus workflow for the ENCODE chip_seq pipeline.')

    ap.add_argument('-i', '--in_files',
                    help='Input files.',
                    nargs='+',
                    required=False)

    ap.add_argument('-p', '--project_name',
                    help='DNAnexus project name',
                    required=True)

    return ap.parse_args()

def upload_files(fl, project, folder):
    files = []
    cmd = 'ua -p ' + project.get_id() + ' -f ' + folder + ' ' + ' '.join(fl)
    print cmd
    output = subprocess.check_output(cmd, shell=True)
    for f in output.split('\n'):
        files += [dxpy.dxlink(f)]

    return files

def get_files(project):
    files = dxpy.find_data_objects(classname='file', name='*.fastq*',
                                   name_mode='glob', project=project.get_id(),
                                   folder=INPUT_FASTQS_FOLDER, return_handler=True)
    files = [dxpy.dxlink(f.get_id()) for f in files]

    return files

def split_reads(fastq_files):
    left_reads = []
    right_reads = []

    for f in fastq_files:
        fn = dxpy.describe(f)['name']
        prefix = re.findall('([\w-]+)', fn)[0]
        if prefix[-1] == '1':
            left_reads += [f]
        elif prefix[-1] == '2':
            right_reads += [f]
        else:
            print 'Can not determine if file is left or right read.  Skipping {0}.'.format(fn)

    return left_reads, right_reads

def populate_workflow(wf, fastq_files, project_name):
    '''This function will populate the workflow for the ChIP-Seq Pipeline.'''
    ref_genome = find_reference_file_by_name(REFERENCE_DATA['reference_genome'])
    ref_genome_index = find_reference_file_by_name(REFERENCE_DATA['reference_genome_bwa_index'])
    (left_reads, right_reads) = split_reads(fastq_files)
    assert len(left_reads) == len(right_reads)

    bwa_input = {'fastq_1': left_reads,
                 'fastq_2': right_reads,
                 'ref_genome': ref_genome,
                 'ref_genome_index': ref_genome_index}
    stage_id = wf.add_stage(find_applet_by_name('bwa'), stage_input=bwa_input, folder='/BWA')
    bwa_output = dxpy.dxlink({'stage': stage_id, 'outputField': 'sorted_bam'})

    pseudo_replicate_input = {}
    stage_id = wf.add_stage()
    pseudo_replicate_output = dxpy.dxlink({'stage': stage_id, 'outputField': 'FIXME'})

    peak_calling_input = {}
    stage_id = wf.add_stage(find_applet_by_name('peak_calling'), stage_input=peak_calling_input, folder='/peak_calling')
    peak_calling_output = dxpy.dxlink({'stage': stage_id, 'outputField': 'FIXME'})

    idr_input = {}
    stage_id = wf.add_stage(find_applet_by_name('idr'), stage_input=idr_input, folder='/idr')
    idr_output = dxpy.dxlink({'stage': stage_id, 'outputField': 'FIXME'})

if __name__ == '__main__':
    args = get_args()

    project = get_project(args.project_name)
    print 'Project: ' + project.describe()['name']

    if args.in_files is not None:
        project.new_folder('/input_fastqs', True)
        fastq_files = upload_files(args.in_files, project, INPUT_FASTQS_FOLDER)
    else:
        fastq_files = get_files(project)

    # Now create a new workflow
    wf = dxpy.new_dxworkflow(title='dx_chip_seq',
                             name='ENCODE ChIP-Seq 2.0',
                             description='The ENCODE ChIP-Seq Pipeline 2.0',
                             project=project.get_id())
    populate_workflow(wf, fastq_files, project.describe()['name'])

