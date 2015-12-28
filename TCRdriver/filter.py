'''methods to filter common contaminants (phix, ribosomal, etc)'''
from __future__ import print_function, division, absolute_import
import sys, os, re, subprocess

BOWTIE2_PHIX_REF = '/Users/mfitzgib/data/Reference/Bowtie2Index/phix'
BOWTIE2_HUMAN_RIBO_REF = '/Users/mfitzgib/data/Reference/Bowtie2Index/human_ribo'

def get_software_version(cmd, version_tag):
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        ret = proc.communicate()
    except OSError:
        return None

    version = None
    for r in ret:
        for line in r.splitlines():
            idx = line.find(version_tag)
            if idx < 0: continue
            version = line[idx + len(version_tag):]
            break

    if version is None:
        return None
    version = version.split()[0]
    return [int(x) for x in version.split('.')]

def get_bowtie_version():
    return get_software_version(['bowtie2', '--version'], ' version ')

def get_pandaseq_version():
    return get_software_version(['pandaseq', '-h'], 'pandaseq ')

def bowtie2_parse_output(ret):
    fields = ret.splitlines()
    total = fields[1]
    unaligned = fields[2]
    assert 'aligned' in unaligned and unaligned.endswith(' 0 times'), unaligned
    total = int(total.split()[0])
    unaligned = int(unaligned.split()[0])
    return total - unaligned, total

def bowtie2_filter(outbase, ref, r1, r2=None):
    outfq = outbase + "_R1.fq.gz"
    cmd = ['bowtie2', '-x', ref]
    if r2 is not None:
        cmd  += ['-1', r1, '-2', r2, '--un-conc-gz', outfq]
    else:
        cmd += ['-U', r1, '--un-gz', outfq]
    cmd += ['-S', '/dev/null']
    proc = subprocess.Popen(cmd, stderr=subprocess.PIPE)
    ret = proc.communicate()
    filtered, total = bowtie2_parse_output(ret[1])

    # rename bowtie2 default paired read output
    if r2 is not None:
        tmp1 = outbase + '_R1.fq.1.gz'
        tmp2 = outbase + '_R1.fq.2.gz'
        out1 = outbase + '_R1.fq.gz'
        out2 = outbase + '_R2.fq.gz'
        os.rename(tmp1, out1)
        os.rename(tmp2, out2)
        outfiles = (out1, out2)
    else:
        outfiles = (outfq, )
    
    return outfiles, filtered, total

def bowtie2_filter_phix(outbase, r1, r2=None):
    return bowtie2_filter(outbase, BOWTIE2_PHIX_REF, r1, r2)
