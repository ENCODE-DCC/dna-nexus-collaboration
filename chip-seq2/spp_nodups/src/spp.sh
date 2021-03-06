#!/bin/bash
# spp 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

main() {
    set -e -x 
    #apt-cache showpkg r-base
    sudo apt-get -f -y --force-yes install r-base-core=2.15.3-1precise0precise1

    echo "Value of input_bam: '$input_bam'"
    echo "Value of control_bam: '$control_bam'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    input_name=`dx describe --name "$input_bam"`
    dx download "$input_bam" -o "$input_name"
    control_name=`dx describe --name "$control_bam"`
    dx download "$control_bam" -o "$control_name"

    pwd
    R --version
    ls /usr/bin/R*
    R CMD INSTALL spp_1.10.1.tar.gz

    # Fill in your application code here.
    #
    # To report any recognized errors in the correct format in
    # $HOME/job_error.json and exit this script, you can use the
    # dx-jobutil-report-error utility as follows:
    #
    #   dx-jobutil-report-error "My error message"
    #
    # Note however that this entire bash script is executed with -e
    # when running in the cloud, so any line which returns a nonzero
    # exit code will prematurely exit the script; if no error was
    # reported in the job_error.json file, then the failure reason
    # will be AppInternalError with a generic error message.

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    Rscript run_spp_nodups.R -c="$input_name" -p=`nproc` -i="$control_name" -npeak=300000 -savr -savp -savd -rf
    
    ls
    
    pdf="${input_name%.*}.pdf"
    rdata="${input_name%.*}.Rdata"
    #peaks=`ls "${input_name%.bam}*.gz"`
    peaks=`ls *.regionPeak.gz`

    file_id=`dx upload "$pdf" --brief`
    dx-jobutil-add-output pdf "$file_id"
    file_id=`dx upload "$rdata" --brief`
    dx-jobutil-add-output rdata "$file_id"
    file_id=`dx upload "$peaks" --brief`
    dx-jobutil-add-output peaks "$file_id"
}
