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
import time
import random

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

# load Eurie's manifest tsv as [Map(path,size,md5),...]
def load_file_list():
    with gzip.open('/ENCODE-SDSC-snapshot-20140505.tsv.gz') as tsvin:
        tsvin = csv.reader(tsvin, delimiter='\t')
        tsvin.next() # skip header
        return [Map(path=row[0], size=int(row[1]), md5=row[2]) for row in tsvin]

# load Jim's patches to Eurie's manifest as [{edwName: Map(badMd5,goodMd5)}]
def load_file_list_patch():
    ans = {}
    with open('/diffMd5.lst') as pin:
        pin = csv.reader(pin, delimiter=' ')
        pin.next() # skip header
        for row in pin:
            ans[row[0]] = Map(path=row[0],badMd5=row[1],goodMd5=row[2])
    print 'Loaded', len(ans), 'MD5 corrections from diffMd5.lst'
    return ans

def disk_free_space():
    s = os.statvfs('/')
    return s.f_bavail * s.f_frsize

# perform the transfer of one file (no retry logic)
def transfer_file(f):
    print f.path, "transfer"

    dn, fn = os.path.split(f.path)
    dn = "/" + dn

    GiB = 1024*1048576
    if f.size > GiB:
        # random delay to ameliorate out-of-the-gate race condition in the
        # following disk free space check
        time.sleep(10*random.random())
    while f.size + GiB > disk_free_space():
        print f.path, "awaiting disk space; need", f.size, "free", disk_free_space()
        time.sleep(60)

    # timeout of last resort -- kill transfer that's gone slower than 4GB/hr (~10 megabit/s)
    last_resort_timeout = max(3600,int(900.0*f.size/GiB))
    timeout_cmd = ["timeout", str(last_resort_timeout)]

    # perform aria2c download
    urls = [base + f.path for base in mirror_urls]
    subprocess.check_call(timeout_cmd
                          + ["aria2c","-x","16", "-s","64","-m","8","--retry-wait=10",
                             "--file-allocation=falloc","-o",f.md5,"-q","--check-certificate=false"]
                          + urls)
    
    try:
        # verify MD5
        md5 = subprocess.Popen(["md5sum",f.md5],stdout=subprocess.PIPE).communicate()[0].split()[0]
        if md5 != f.md5:
            raise Exception("expected MD5 " + f.md5 + " got " + md5)
        print f.path, "verified", md5, f.md5
        
        # run DNAnexus upload agent to put the file in the project
        subprocess.check_call(timeout_cmd
                              + ["/ua","-p",dxpy.PROJECT_CONTEXT_ID,"-f",dn,"-n",fn,"-r","8",
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
    print f.path, "begin", f.size, f.md5
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
            print f.path, "removing incomplete file", existing_dxfile.get_id()
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
            if n >= max_tries-1:
                raise
            print f.path, "retry"

    print f.path, "complete", f.size
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
    dstat = subprocess.Popen(['dstat','-cmnd','--freespace','60'])
    try:
        subprocess.check_call("gunzip /ua.gz; chmod +x /ua",shell=True)
        files = load_file_list()
        files_patch = load_file_list_patch()

        # sort the files by MD5 (effectively shuffle them, but such that the
        # order is the same across parallel workers)
        if not smallest:
            files.sort(key=lambda f: f.md5)
        else:
            # for rapid testing: start with the smallest files
            files.sort(key=lambda f: f.size)

        # take a subset for this worker
        files = [files[i] for i in xrange(len(files)) if i%workers == whoami]
        # apply MD5 patches
        for f in files:
            if f.path in files_patch:
                if f.md5 != files_patch[f.path].badMd5:
                    raise Exception(f.path, "bad MD5 patch", f.md5, files_patch[f.path].badMd5, files_patch[f.path].goodMd5)
                f.md5 = files_patch[f.path].goodMd5
                print f.path, "patched MD5 from", files_patch[f.path].badMd5, "to", f.md5
        for f in files:
            if f.path == "2013/4/18/ENCFF001QYA.bam":
                print "For avoidance of doubt, expected MD5 for 2013/4/18/ENCFF001QYA.bam is now", f.md5
        # apply max_files_per_worker
        if max_files_per_worker is not None and len(files) > max_files_per_worker:
            del files[max_files_per_worker:]

        print "will process", len(files), "files totaling", sum([f.size for f in files]), "bytes"

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
def main(workers, max_files_per_worker=None, threads_per_worker=8, worker_launch_delay_seconds=0, smallest=False):
    mkdirs()

    worker_instance_type = "mem2_hdd2_x4"
    if smallest:
        # debugging - run on default instances
        worker_instance_type = None

    # launch workers, each to process a subset of the files
    subjobs = []
    for i in range(workers):
        subjob_input = { "workers": workers, "max_files_per_worker": max_files_per_worker, "whoami": i,
                         "threads_per_worker": threads_per_worker, "smallest": smallest }
        subjobs.append(dxpy.new_dxjob(subjob_input, "process", instance_type=worker_instance_type))
        if worker_launch_delay_seconds > 0 and i < (workers-1):
            # delay launching each worker to smooth out the load on the remote
            # server
            time.sleep(worker_launch_delay_seconds)

    # schedule postprocessing to reduce statistics
    output_fields = ["files_skipped", "files_transferred", "bytes_transferred"]
    postprocess_job = dxpy.new_dxjob(fn_input={k:[subjob.get_output_ref(k) for subjob in subjobs] for k in output_fields},
                                     fn_name="postprocess")

    return {k:postprocess_job.get_output_ref(k) for k in output_fields}

dxpy.run()
