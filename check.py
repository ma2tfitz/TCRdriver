#!/usr/bin/env python

from __future__ import print_function, division, absolute_import
import sys
import os
import argparse
import TCRdriver

r1 = "./data/150504-subset_R1.fastq.gz"
r2 = "./data/150504-subset_R2.fastq.gz"

run_base = TCRdriver.fastq_basename(r1)
if r2 is not None:
    assert TCRdriver.fastq_basename(r2) == run_base
print("Run name:", run_base)
nophix_base = run_base + '-nophix'

nophix_files, phix, total = TCRdriver.filter.bowtie2_filter_phix(nophix_base, r1, r2)
print(total, "total reads")
print(phix/total, "phiX percentage")

merge_reads = False

if merge_reads:
    merged_reads = TCRdriver.merge.pandaseq(run_base, *nophix_files)
    print("merged", merged_reads)
    total, counter = TCRdriver.demux.demux("150504-DavidPilot-barcodes.tsv", merged_reads)
    for sample in counter:
        if sample == "unmatched": continue
        r1 = "{0}_R1.fq.gz".format(sample)
        TCRdriver.tcr.mixcr_rna(sample, r1)
else:
    total, counter = TCRdriver.demux.demux("150504-DavidPilot-barcodes.tsv", *nophix_files)
    for sample in counter:
        if sample == "unmatched": continue
        r1 = "{0}_R1.fq.gz".format(sample)
        r2 = "{0}_R2.fq.gz".format(sample)
        TCRdriver.tcr.mixcr_rna(sample, r1, r2)
