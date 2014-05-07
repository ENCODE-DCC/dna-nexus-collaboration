#!/usr/bin/env python
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import os
import dxpy
import gzip
import csv
import subprocess
import threading
from multiprocessing.dummy import Pool as ThreadPool

# URL prefixes of data mirrors. aria2c can download using parallel connections
# to multiple mirrors.
mirror_urls = [
    "http://encodedcc.sdsc.edu/warehouse/",
    "http://encode-01.sdsc.edu/warehouse/"
]

# utility class
class Map(dict):
    def __init__(self, **kwargs):
        super(Map, self).__init__(**kwargs)
        self.__dict__ = self

# load Eurie's manifest tsv as [Map(filepath,size,md5),...]
def load_file_list():
    with gzip.open('/ENCODE-SDSC-snapshot-20140505.tsv.gz') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        tsvin.next() # skip header
        return [Map(path=row[0], size=int(row[1]), md5=row[2]) for row in tsvin]

# perform the transfer of one file (no retry logic)
def transfer_file(f):
    print f.path, "transfer"

    dn, fn = os.path.split(f.path)
    dn = "/" + dn

    # perform aria2c download
    urls = [base + f.path for base in mirror_urls]
    subprocess.check_call(["aria2c","-x","16", "-s","64","-m","8","--retry-wait=10",
                           "-o",f.md5,"-q","--check-certificate=false"] + urls)
    
    try:
        # verify MD5
        md5 = subprocess.Popen(["md5sum",f.md5],stdout=subprocess.PIPE).communicate()[0].split()[0]
        if md5 != f.md5:
            raise Exception("expected MD5 " + f.md5 + " got " + md5)
        print f.path, "verified", md5, f.md5
        
        # run DNAnexus upload agent to put the file in the project
        subprocess.check_call(["/ua","-p",dxpy.PROJECT_CONTEXT_ID,"-f",dn,"-n",fn,"-r","8",
                               "--do-not-compress","--do-not-resume",f.md5])
        
        # set md5 property on the new file object
        dxfiles = list(dxpy.search.find_data_objects(project=dxpy.PROJECT_CONTEXT_ID,classname="file",
                                                     folder=dn,name=fn,return_handler=True))
        if len(dxfiles) != 1:
            raise Exception("could not uniquely locate file just uploaded")
        dxfiles[0].set_properties({"md5": f.md5})
    finally:
        try:
            # delete scratch copy
            os.remove(f.md5)
        except:
            pass

# "process" one file idempotently -- check if it's already present in the
# project, and if not, perform the transfer with retry logic
def process_file(f):
    print f.path, "begin"
    dn, fn = os.path.split(f.path)
    dn = "/" + dn

    # check if a file with this path & MD5 already exists in the project
    existing_dxfiles = list(dxpy.search.find_data_objects(project=dxpy.PROJECT_CONTEXT_ID,classname="file",
                                                          folder=dn,name=fn,return_handler=True))
    if len(existing_dxfiles) > 1:
        raise Exception(f.path + " found multiple existing file objects! manually remove them and try again")
    elif len(existing_dxfiles) == 1:
        existing_dxfile = existing_dxfiles[0]
        if existing_dxfile.state == "open":
            print f.path, " removing incomplete file", existing_dxfile.get_id()
            existing_dxfile.remove()
        else:
            existing_props = existing_dxfile.get_properties()
            if "md5" not in existing_props or existing_props["md5"] != f.md5:
                raise Exception(f.path + " MD5 mismatch with existing file object! manually remove and try again")
            print f.path, "skip", f.md5, existing_props["md5"], existing_dxfile.get_id()
            return 0

    # perform transfer, with retry logic
    max_tries = 8
    for n in xrange(max_tries):
        try:
            transfer_file(f)
            break
        except:
            print f.path, sys.exc_info()
            if n == max_tries-1:
                raise
            print f.path, "retry"

    print f.path, "complete"
    return f.size

# ensure existence of all necessary folders in the project
def mkdirs():
    files = load_file_list()
    dirs = set([os.path.dirname(f.path) for f in files])

    proj = dxpy.DXProject(dxpy.PROJECT_CONTEXT_ID)
    for dn in dirs:
        proj.new_folder("/" + dn, parents=True)

@dxpy.entry_point("process")
def process(workers, max_files_per_worker, whoami, threads_per_worker, smallest):
    dstat = subprocess.Popen(['dstat','-cmdn','60'])
    try:
        subprocess.check_call("gunzip /ua.gz; chmod +x /ua",shell=True)
        files = load_file_list()

        # sort the files by MD5 (effectively shuffle them, but such that the
        # order is the same across parallel workers)
        if not smallest:
            files.sort(key=lambda f: f.md5)
        else:
            # for rapid testing: start with the smallest files
            files.sort(key=lambda f: f.size)

        # take a subset for this worker
        files = [files[i] for i in xrange(len(files)) if i%workers == whoami]
        if len(files) > max_files_per_worker:
            del files[max_files_per_worker:]

        print "will process", len(files), "files"

        # launch a thread pool to process them
        pool = ThreadPool(threads_per_worker)
        results = pool.map(process_file, files)

        # compute stats...assumes there are no legit zero-length files
        files_skipped = len([x for x in results if x == 0])
        files_transferred = len([x for x in results if x > 0])
        bytes_transferred = sum(results)
        print "all done, transferred", files_transferred, "files and", bytes_transferred, "bytes; skipped", files_skipped, "files"

        return { "files_skipped": files_skipped, "files_transferred": files_transferred,
                 "bytes_transferred": bytes_transferred }
    finally:
        dstat.kill()

@dxpy.entry_point("postprocess")
def postprocess(files_skipped, files_transferred, bytes_transferred):
    total_files_skipped = sum(files_skipped)
    total_files_transferred = sum(files_transferred)
    total_bytes_transferred = sum(bytes_transferred)
    return { "files_skipped": total_files_skipped, "files_transferred": total_files_transferred,
             "bytes_transferred": total_bytes_transferred }

@dxpy.entry_point("main")
def main(workers, max_files_per_worker=None, threads_per_worker=8, smallest=False):
    mkdirs()

    # launch workers, each to process a subset of the files
    subjobs = []
    for i in range(workers):
        subjob_input = { "workers": workers, "max_files_per_worker": max_files_per_worker, "whoami": i,
                         "threads_per_worker": threads_per_worker, "smallest": smallest }
        subjobs.append(dxpy.new_dxjob(subjob_input, "process"))

    # schedule postprocessing to reduce statistics
    output_fields = ["files_skipped", "files_transferred", "bytes_transferred"]
    postprocess_job = dxpy.new_dxjob(fn_input={k:[subjob.get_output_ref(k) for subjob in subjobs] for k in output_fields},
                                     fn_name="postprocess")

    return {k:postprocess_job.get_output_ref(k) for k in output_fields}

dxpy.run()
