"""methods to process TCR fastq files"""
from __future__ import print_function, division, absolute_import
import sys
import os
import re
import subprocess

def mixcr_rna_single(outbase, r1):
    
    if outbase.endswith("TRA"):
        locus = "TRA"
    elif outbase.endswith("TRB"):
        locus = "TRB"
    else:
        raise Except

    alignment_filename = outbase + "_mixcr.vdjca"
    clone_filename = outbase + "_mixcr.clns"
    export_filename = outbase + "_mixcr.tsv"

    # not using "-p rna-seq"  for now
    cmd = ["mixcr", "align", "--loci", locus, "--species", "hs", r1, alignment_filename]
    print("running:", *cmd)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    print(outbase, ret)

    cmd = ["mixcr", "assemble", alignment_filename, clone_filename]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    print(outbase, ret)

    cmd = ["mixcr", "exportClones", "-s", clone_filename, export_filename]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    print(outbase, ret)

def mixcr_rna_paired(outbase, r1, r2):
    if outbase.endswith("TRA"):
        locus = "TRA"
    elif outbase.endswith("TRB"):
        locus = "TRB"
    else:
        raise Except

    alignment_filename = outbase + "_mixcr.vdjca"
    clone_filename = outbase + "_mixcr.clns"
    export_filename = outbase + "_mixcr.tsv"

    # not using "-p rna-seq"  for now
    cmd = ["mixcr", "align", "--loci", locus, "--species", "hs", r1, r2, alignment_filename]
    print("running:", *cmd)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    print(outbase, ret)

    cmd = ["mixcr", "assemble", alignment_filename, clone_filename]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    print(outbase, ret)

    cmd = ["mixcr", "exportClones", "-s", clone_filename, export_filename]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    print(outbase, ret)

def mixcr_rna(outbase, r1, r2=None):
    
    if outbase.endswith("TRA"):
        locus = "TRA"
    elif outbase.endswith("TRB"):
        locus = "TRB"
    else:
        raise Except

    alignment_filename = outbase + "_mixcr.vdjca"
    clone_filename = outbase + "_mixcr.clns"
    export_filename = outbase + "_mixcr.tsv"

    # not using "-p rna-seq"  for now
    cmd = ["mixcr", "align", "-f", "--loci", locus, "--species", "hs"]
    if r2 is None:
        cmd += [r1, ]
    else:
        cmd += [r1, r2, ]
    cmd += [alignment_filename, ]
    print("running:", *cmd)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    print(outbase, ret)

    cmd = ["mixcr", "assemble", "-f", alignment_filename, clone_filename]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    print(outbase, ret)

    cmd = ["mixcr", "exportClones", "-s", clone_filename, export_filename]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    ret = proc.communicate()
    print(outbase, ret)

def mitcr_rna_single(outbase, r1):
    assert False, "mitcr_rna_single not yet implemented"
