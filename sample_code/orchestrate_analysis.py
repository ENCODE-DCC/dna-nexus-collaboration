#!/usr/bin/env python

import sys
import json
import dxpy

FACTORY_PROJECT = 'project-BJ4517Q0PZbjk7K76v1004FP'
COMPLETED_ANALYSIS = 'analysis-BJ465j00PZbX8JV9K59003qf'

def main(token):
    # Configure dxpy authentication
    dxpy.set_security_context({'auth_token_type': 'Bearer', 'auth_token': token})

    # Resolve FACTORY_PROJECT by ID
    proj = dxpy.DXProject(FACTORY_PROJECT)
    print 'Resolved project:', proj.describe()['name'], proj.get_id()

    # Set FACTORY_PROJECT as the workspace for subsequent operations
    # (sort of like the current working directory)
    dxpy.set_workspace_id(FACTORY_PROJECT)

    # Resolve the workflow by name. (Could also store ID like the project)
    wf = list(dxpy.search.find_data_objects(classname="workflow", name="RNA-seq pipeline",
                                            return_handler=True))[0]
    print 'Resolved workflow:', wf.describe()['name'], wf.get_id()

    # TODO: Stage the inputs. Here we find them in the IN folder
    left_reads = list(dxpy.search.find_data_objects(classname="file", name="ENCFF001JPX.1k.fastq.gz",
                                                    folder="/IN", return_handler=True))[0]
    print 'Resolved left reads:', left_reads.describe()['name'], left_reads.get_id()
    right_reads = list(dxpy.search.find_data_objects(classname="file", name="ENCFF001JQB.1k.fastq.gz",
                                                     folder="/IN", return_handler=True))[0]
    print 'Resolved right reads:', right_reads.describe()['name'], right_reads.get_id()

    # Launch the workflow
    analysis = wf.run({'0.fastqs': [dxpy.dxlink(left_reads.get_id())],
                       '0.fastq_pairs': [dxpy.dxlink(right_reads.get_id())]})
    print 'Launched analysis:', analysis.get_id()
    print 'Analysis state:', analysis.describe()['state']

    # TODO: Poll for (or come back when) analysis state 'done' or 'failed'.
    # Handle any failures.

    # Cooking-show-style substitution with completed analysis
    analysis = dxpy.DXAnalysis(COMPLETED_ANALYSIS)
    print 'Analysis state:', analysis.describe()['state']

    # Enumerate outputs
    print 'Analysis outputs:'
    for one_output_name, one_output_link in analysis.describe()['output'].iteritems():
        one_output = dxpy.get_handler(one_output_link) # one_output : dxpy.DXFile
        one_file_name = one_output.describe()['name']
        one_file_url, _ = one_output.get_download_url(preauthenticated=True, filename=one_file_name)
        print one_file_name, one_file_url

if len(sys.argv) < 2:
    print 'Usage: orchestrate_analysis.py <auth_token>'
    sys.exit(1)

main(sys.argv[1])
