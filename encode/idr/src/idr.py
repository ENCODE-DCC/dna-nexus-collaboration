#!/usr/bin/env python
# idr 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# See https://wiki.dnanexus.com/Developer-Portal for documentation and
# tutorials on how to modify this file.
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import os
import itertools
import subprocess
import multiprocessing

import dxpy

def download_and_gunzip_file(input_file, skip_decompress=False):
    input_file = dxpy.DXFile(input_file)
    input_filename = input_file.describe()['name']
    ofn = input_filename

    cmd = 'dx download ' + input_file.get_id() + ' -o '
    if input_filename.endswith('.tar.gz') and not skip_decompress:
        cmd += '- | tar -zxvf - > '
        ofn = ofn.replace('.tar.gz', '')
    elif (os.path.splitext(input_filename)[-1] == '.gz') and not skip_decompress:
        cmd += '- | gunzip > '
        ofn = os.path.splitext(ofn)[0]
    cmd += ofn
    print cmd
    subprocess.check_call(cmd, shell=True)

    return ofn

def download_dx_files(dx_files, skip_decompress=False):
    results = []
    for dx_file in dx_files:
        results += [download_and_gunzip_file(dx_file, skip_decompress)]
    return results

def setup_idr_batch_consistency_analysis(job_inputs):
    # Copy stuff from the /opt/idrCode directory to satisfy the hard-coded
    # paths in batch-consistency-analysis.r.
    print 'genome_table_filename: {0}'.format(job_inputs['genome_table_filename'])

    cmd = ['cp', '/opt/idrCode/functions-all-clayton-12-13.r', '.']
    print cmd
    subprocess.check_call(cmd)
    cmd = ['cp', '/opt/idrCode/genome_tables/' + job_inputs['genome_table_filename'], 'genome_table.txt']
    print cmd
    subprocess.check_call(cmd)

def idr_batch_consistency_analysis(job_inputs):
    """Runs IDR batch consistency analysis on two files containing peaks.
    Outputs several files produced by the IDR R script.
    """
    print 'peaks_file1: {0}'.format(job_inputs['peaks_file1'])
    print 'peaks_file2: {0}'.format(job_inputs['peaks_file2'])
    print 'output_prefix: {0}'.format(job_inputs['output_prefix'])
    print 'ranking_measure: {0}'.format(job_inputs['ranking_measure'])

    # Run IDR
    cmd = ['Rscript', '/opt/idrCode/batch-consistency-analysis.r', job_inputs['peaks_file1'], job_inputs['peaks_file2'],
           '-1', job_inputs['output_prefix'], '0', 'F', job_inputs['ranking_measure']]
    print cmd
    subprocess.check_call(cmd)

    cmd = 'tar zcvf {0}.tar.gz {0}-em.sav {0}-uri.sav {0}-Rout.txt {0}-npeaks-aboveIDR.txt {0}-overlapped-peaks.txt'.format(job_inputs['output_prefix'])
    print cmd
    subprocess.check_call(cmd, shell=True)

    return '{0}.tar.gz'.format(job_inputs['output_prefix'])

def get_peak_counts(fn, threshold):
    cmd = "awk '$11 <= " + str(threshold) + " {print $0}' " + fn + " | wc -l "
    print cmd
    peak_counts = subprocess.check_output(cmd, shell=True)

    return int(peak_counts)

def get_thresholds(replicate_idr_prefixes, pseudo_replicate_idr_prefixes, pooled_pseudo_replicate_idr_prefix, replicate_peaks_threshold, pseudo_replicate_peaks_threshold, pooled_pseudo_replicate_peaks_threshold):
    replicate_thresholds = []
    for replicate in replicate_idr_prefixes:
        replicate_thresholds += [get_peak_counts(replicate + '-overlapped-peaks.txt', replicate_peaks_threshold)]

    pseudo_replicate_thresholds = []
    for pseudo_replicate in pseudo_replicate_idr_prefixes:
        pseudo_replicate_thresholds += [get_peak_counts(pseudo_replicate + '-overlapped-peaks.txt', pseudo_replicate_peaks_threshold)]

    pooled_pseudo_replicate_thresholds = get_peak_counts(pooled_pseudo_replicate_idr_prefix + '-overlapped-peaks.txt', pooled_pseudo_replicate_peaks_threshold)


    num_peaks_each_rep = replicate_thresholds
    num_peaks_each_pseudo_rep = pseudo_replicate_thresholds
    numPeaks_Rep0 = pooled_pseudo_replicate_thresholds

    return num_peaks_each_rep, num_peaks_each_pseudo_rep, numPeaks_Rep0

def create_final_set_of_peak_calls(job_inputs):
    replicate_idr_prefixes = [r.replace('.tar.gz', '') for r in job_inputs['replicate_idr_files']]
    pseudo_replicate_idr_prefixes = [r.replace('.tar.gz', '') for r in job_inputs['pseudo_replicate_idr_files']]
    pooled_pseudo_replicate_idr_prefix = job_inputs['pooled_pseudo_replicate_idr_files'].replace('.tar.gz', '')

    (num_peaks_each_rep, num_peaks_each_pseudo_rep, numPeaks_Rep0) = get_thresholds(replicate_idr_prefixes,
                                                                                    pseudo_replicate_idr_prefixes,
                                                                                    pooled_pseudo_replicate_idr_prefix,
                                                                                    job_inputs['replicate_peaks_threshold'],
                                                                                    job_inputs['pseudo_replicate_peaks_threshold'],
                                                                                    job_inputs['pooled_pseudo_replicate_peaks_threshold'])
    max_numPeaks_Rep = max(num_peaks_each_rep)

    pooled_replicates_peaks_fn = download_and_gunzip_file(job_inputs['pooled_replicate_peaks_file'])
    coi = {'signal.value': 7, 'p.value': 8, 'q.value': 9}[job_inputs['ranking_measure']]
    cmd = 'sort -k{0}nr,{0}nr {1} | head -n {2} | gzip -c > {3}_conservative.regionPeak.gz'.format(coi, pooled_replicates_peaks_fn, max_numPeaks_Rep, job_inputs['output_prefix'])
    print cmd
    subprocess.check_output(cmd, shell=True)

    opt_thresh = max(max_numPeaks_Rep, numPeaks_Rep0)
    cmd = 'sort -k{0}nr,{0}nr {1} | head -n {2} | gzip -c > {3}_optimal.regionPeak.gz'.format(coi, pooled_replicates_peaks_fn, opt_thresh, job_inputs['output_prefix'])
    print cmd
    subprocess.check_output(cmd, shell=True)

    conservative_result = dxpy.upload_local_file('{0}_conservative.regionPeak.gz'.format(job_inputs['output_prefix']))
    optimal_result = dxpy.upload_local_file('{0}_optimal.regionPeak.gz'.format(job_inputs['output_prefix']))

    return {'conservative_peak_calls': dxpy.dxlink(conservative_result),
            'optimal_peak_calls': dxpy.dxlink(optimal_result),
            'num_peaks_each_rep': num_peaks_each_rep,
            'num_peaks_each_pseudo_rep': num_peaks_each_pseudo_rep,
            'num_peaks_pooled_pseudo_rep': numPeaks_Rep0}

@dxpy.entry_point('main')
def main(**job_inputs):
    pool = multiprocessing.Pool()

    setup_idr_batch_consistency_analysis(job_inputs)

    replicate_peaks_filenames = download_dx_files(job_inputs['replicate_peaks_files'])
    pseudo_replicate_peaks_filenames = download_dx_files(job_inputs['pseudo_replicate_peaks_files'])
    pooled_pseudo_replicate_filenames = download_dx_files(job_inputs['pooled_pseudo_replicate_peaks_file'])

    replicate_idr_output = []
    assert len(job_inputs['replicate_peaks_files']) > 1
    for i, j in itertools.combinations(range(len(replicate_peaks_filenames)), 2):
        prefix = job_inputs['output_prefix'] + '_replicate_{0}_vs_{1}'.format(i, j)
        idr_batch_consistency_analysis_input = {'peaks_file1': replicate_peaks_filenames[i],
                                                'peaks_file2': replicate_peaks_filenames[j],
                                                'output_prefix': prefix,
                                                'ranking_measure': job_inputs['ranking_measure'],
                                                'genome_table_filename': job_inputs['genome_table_filename']}
        replicate_idr_output += [pool.apply_async(idr_batch_consistency_analysis, (idr_batch_consistency_analysis_input, ))]


    pseudo_replicate_output = []
    assert (len(pseudo_replicate_peaks_filenames) % 2) == 0
    for i in xrange(0, len(pseudo_replicate_peaks_filenames), 2):
        prefix = job_inputs['output_prefix'] + '_pseudo_replicate_{0}'.format(i)
        idr_batch_consistency_analysis_input = {'peaks_file1': pseudo_replicate_peaks_filenames[i],
                                                'peaks_file2': pseudo_replicate_peaks_filenames[i+1],
                                                'output_prefix': prefix,
                                                'ranking_measure': job_inputs['ranking_measure'],
                                                'genome_table_filename': job_inputs['genome_table_filename']}
        pseudo_replicate_output += [pool.apply_async(idr_batch_consistency_analysis, (idr_batch_consistency_analysis_input, ))]

    assert len(pooled_pseudo_replicate_filenames) == 2
    prefix = job_inputs['output_prefix'] + '_pooled_pseudo_replicate'
    idr_batch_consistency_analysis_input = {'peaks_file1': pooled_pseudo_replicate_filenames[0],
                                            'peaks_file2': pooled_pseudo_replicate_filenames[1],
                                            'output_prefix': prefix,
                                            'ranking_measure': job_inputs['ranking_measure'],
                                            'genome_table_filename': job_inputs['genome_table_filename']}
    pooled_pseudo_replicate_output = pool.apply_async(idr_batch_consistency_analysis, (idr_batch_consistency_analysis_input, ))

    pool.close()
    pool.join()
    replicate_idr_output = [r.get(timeout=10) for r in replicate_idr_output]
    pseudo_replicate_output = [r.get(timeout=10) for r in pseudo_replicate_output]
    pooled_pseudo_replicate_output = pooled_pseudo_replicate_output.get(timeout=10)

    create_final_set_of_peak_calls_input = job_inputs.copy()
    create_final_set_of_peak_calls_input['replicate_idr_files'] = replicate_idr_output
    create_final_set_of_peak_calls_input['pseudo_replicate_idr_files'] = pseudo_replicate_output
    create_final_set_of_peak_calls_input['pooled_pseudo_replicate_idr_files'] = pooled_pseudo_replicate_output
    create_final_set_of_peak_calls_output = create_final_set_of_peak_calls(create_final_set_of_peak_calls_input)

    output = {'conservative_peak_calls': create_final_set_of_peak_calls_output['conservative_peak_calls'],
              'optimal_peak_calls': create_final_set_of_peak_calls_output['optimal_peak_calls'],
              'num_peaks_each_rep': create_final_set_of_peak_calls_output['num_peaks_each_rep'],
              'num_peaks_each_pseudo_rep': create_final_set_of_peak_calls_output['num_peaks_each_pseudo_rep'],
              'num_peaks_pooled_pseudo_rep': create_final_set_of_peak_calls_output['num_peaks_pooled_pseudo_rep']}

    return output

dxpy.run()
