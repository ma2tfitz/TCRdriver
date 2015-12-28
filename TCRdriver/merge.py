"""methods to merge overlapping paired reads"""
from __future__ import print_function, division, absolute_import
import sys
import os
import re
import subprocess

def pandaseq(outbase, r1, r2):
    unpaired_filename = outbase + "-unpaired.fq"
    log_filename = outbase + "-pandaseq.log"
    output_filename = outbase + "-merge.fq"
    
    cmd = ["pandaseq", "-f", r1, "-r", r2, "-B", "-F", "-T", "8", "-t", ".5", "-U", unpaired_filename, "-d", "rbFkms", "-g", log_filename, "-w", output_filename]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    
    return output_filename
