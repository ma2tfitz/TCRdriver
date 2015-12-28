from __future__ import print_function, division, absolute_import
import sys
import os
import re
import subprocess
import TCRdriver.filter
import TCRdriver.demux
import TCRdriver.merge
import TCRdriver.tcr

version = TCRdriver.filter.get_bowtie_version()
if version is None:
    raise Exception('required software bowtie2 not found')
if not (version[0] >= 2 and version[1] >= 2):
    raise Exception('bowtie2 minimum supported version is 2.2.0')


version = TCRdriver.filter.get_pandaseq_version()
if version is None:
    raise Exception('required software pandaseq not found')
if not (version[0] >= 2 and version[1] >= 7):
    raise Exception('pandaseq minimum supported version is 2.7')

def fastq_basename(fastq):
    f = os.path.basename(fastq)
    for s in ('.gz', '.fq', '.fastq', '_001', '_R1', '_R2'):
        if f.endswith(s):
            f = f[:-len(s)]
    return f



