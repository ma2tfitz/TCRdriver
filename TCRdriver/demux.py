'''simple demux for multi-part barcodes'''
from __future__ import print_function, division, absolute_import

import sys
import os
import re
import subprocess
import gzip

from collections import namedtuple, Counter
from string import maketrans

def revcomp(seq, t=maketrans("ACGTacgt", "TGCAtgca")):
    """very simple reverse complement"""
    return seq.translate(t)[::-1]

def fastqreader(filename):
    infile = gzip.open(filename) if filename.endswith(".gz") else open(filename)
    fields = []
    for line in infile:
        fields.append(line.rstrip())
        if len(fields) >= 4:
            yield fields[0][1:], fields[1], fields[3]
            fields = []
    infile.close()

class SampleName(object):
    def __init__(self, sample, chain, replicate):
        self.sample = sample
        self.chain = chain
        self.replicate = int(replicate)
        self.name = "%s_%s_%s" % (self.sample, self.replicate, self.chain)

class BarcodeDefinition(object):
    def __init__(self, start, end, name):
        assert (start > 0 and end > start) or (start < 0 and end > start)
        self.start = start
        self.end = end
        self.name = name
    def __repr__(self):
        return "BarcodeDefinition(start={0}, end={1}, name={2})".format(self.start, self.end, self.name)

class BarcodeAssignment(object):
    """
    associate a sample name (sample, TCR chain, replicate number)
    with a tuple of sequences expected at positions from the
    barcode_definitions.
    """
    def __init__(self, sample, chain, replicate, sequences):
        self.sample = sample
        self.chain = chain
        self.replicate = replicate
        self.name = "{0}-{1}-{2}".format(sample, replicate, chain)
        self.sequences = sequences

    def __repr__(self):
        return "BarcodeAssignment(name={0}, sequences={1})".format(self.name, self.sequences)
        
class BarcodeCollection(object):
    """
    Collection of barcode definitions (expected locations in sequence)
    and assignments of specific sequences to samples
    """
    def __init__(self):
        self.barcode_definitions = []
        self.barcode_assignments = []

    def _parse_tsv_defline(self, defline):
        """
        Parse the defline in TSV config file
        """
        fields = defline.split("\t")
        assert tuple(fields[:3]) == ("Sample", "Chain", "Replicate"), ",".join(fields[:3])
        for b in fields[3:]:
            if b == "Description":
                break
            i1 = b.find("_")
            i2 = b.find(":", i1)
            assert i1 >= 0 and i2 > i1, fields
            name = b[:i1]
            start = int(b[i1+1:i2])
            end = int(b[i2+1:])
            self.barcode_definitions.append(BarcodeDefinition(start, end, name))

    def parse_tsv(self, config_filename):
        """
        Parse a simple TSV file defining barcodes
        """
        for line in open(config_filename):
            line = line.rstrip()
            if line == "" or line[0] == "#": continue
            if not self.barcode_definitions:
                self._parse_tsv_defline(line)
                n = 3 + len(self.barcode_definitions)
            else:
                fields = line.split("\t")
                sequences = tuple(fields[3:n])
                # TODO check that sequences contain only ACGT here
                assert len(sequences) == len(self.barcode_definitions)
                self.barcode_assignments.append(BarcodeAssignment(fields[0], fields[1], int(fields[2]), sequences))

def _hamming(seq1, seq2):
    """
    Simple hamming distance; assumes sequence lengths are the same.
    Alternatives (eg. Levenshtein.hamming) might save ~5% on total
    runtime.
    """
    if seq1 == seq2:
        return 0
    return sum([s1 != s2 for s1, s2 in zip(seq1, seq2)])

def closest_barcode(seq, barcode_lookup, max_mismatch=1, dist_func=_hamming):
    """
    Find closest single barcode within some allowable number of
    mismatches.  Allowing more than one mismatch or using more complex
    metric (eg. Levenshtein dist) seem to cause more trouble for a
    modest increase in matches.  """
    if seq in barcode_lookup:
        return 0, seq

    if max_mismatch == 0:
        return None, None

    hits = set()
    current_mismatch = 9999
    for candidate in barcode_lookup:
        m = 0
        for s, c in zip(seq, candidate):
            m += dist_func(s, c)
            #if m > max_mismatch: break
        if m <= max_mismatch:
            if m < current_mismatch:
                hits = set()
                current_mismatch = m
            hits.add((m, candidate))

    if len(hits) == 1:
        return hits.pop()

    return None, None

def demux_single(barcode_collection, r1, remove_barcodes=True):
    total = 0
    counter = Counter()
    barcode_lookup = {}

    unmatched_key = "unmatched"

    # setup a hash of output files for writing
    # TODO: this will NOT scale much beyond two plates (96 * 2 * 2 descriptors)
    outfiles = {}
    for b in barcode_collection.barcode_assignments:
        barcode_lookup[b.sequences] = b
        o1 = gzip.open("{0}_R1.fq.gz".format(b.name), "w")
        outfiles[b.name] = o1
    o1 = gzip.open("unmatched_R1.fq.gz", "w")
    outfiles[unmatched_key] = o1

    s1offset = max(0, max([b.end for b in barcode_collection.barcode_definitions]))
    s2offset = min(0, min([b.start for b in barcode_collection.barcode_definitions]))
    if s2offset == 0:
        s2offset = None # for slicing

    for b in barcode_collection.barcode_assignments:
        s = tuple([x if d.start >= 0 else revcomp(x) for d,x in zip(barcode_collection.barcode_definitions, b.sequences)])
        barcode_lookup[s] = b

    for hdr, seq, qual in fastqreader(r1):
        total += 1
        s = tuple([seq[b.start:b.end] for b in barcode_collection.barcode_definitions])
        mismatch, key = closest_barcode(s, barcode_lookup, max_mismatch=1)
        if key is None:
            counter[unmatched_key] += 1
            print("@{0}\n{1}\n+\n{2}".format(hdr, seq, qual), file=outfiles[unmatched_key])
        else:
            b = barcode_lookup[key]
            counter[b.name] += 1
            if remove_barcodes:
                seq = seq[s1offset:s2offset]
                qual = qual[s1offset:s2offset]
            if mismatch > 0:
                hdr = "{0} mismatch={1}".format(hdr, mismatch)
            print("@{0}\n{1}\n+\n{2}".format(hdr, seq, qual), file=outfiles[b.name])

    for key in outfiles:
        outfiles[key].close()

    return total, counter

def demux_pair(barcode_collection, r1, r2, remove_barcodes=True):
    total = 0
    counter = Counter()
    barcode_lookup = {}
    
    unmatched_key = "unmatched"

    # setup a hash of output files for writing
    # TODO: this will NOT scale much beyond one plate (96 * 2 *2 descriptors)
    outfiles = {}
    for b in barcode_collection.barcode_assignments:
        barcode_lookup[b.sequences] = b
        o1 = gzip.open("{0}_R1.fq.gz".format(b.name), "w")
        o2 = gzip.open("{0}_R2.fq.gz".format(b.name), "w")
        outfiles[b.name] = (o1, o2)
    o1 = gzip.open("unmatched_R1.fq.gz", "w")
    o2 = gzip.open("unmatched_R2.fq.gz", "w")
    outfiles[unmatched_key] = (o1, o2)

    s1offset = max(0, max([b.end for b in barcode_collection.barcode_definitions]))
    s2offset = max(0, -min([b.start for b in barcode_collection.barcode_definitions]))

    f1 = fastqreader(r1)
    f2 = fastqreader(r2)

    while True:
        try:
            h1, s1, q1 = f1.next()
            h2, s2, q2 = f2.next()
        except StopIteration:
            break
        
        total += 1

        s = tuple([s1[b.start:b.end] if b.start >= 0 else s2[-b.end:-b.start] for b in barcode_collection.barcode_definitions])
        mismatch, key = closest_barcode(s, barcode_lookup, max_mismatch=1)
        if key is None:
            counter[unmatched_key] += 1
            print("@{0}\n{1}\n+\n{2}".format(h1, s1, q1), file=outfiles[unmatched_key][0])
            print("@{0}\n{1}\n+\n{2}".format(h2, s2, q2), file=outfiles[unmatched_key][1])
        else:
            b = barcode_lookup[key]
            counter[b.name] += 1
            if remove_barcodes:
                s1 = s1[s1offset:]
                q1 = q1[s1offset:]
                s2 = s2[s2offset:]
                q2 = q2[s2offset:]
            if mismatch > 0:
                h1 = "{0} mismatch={1}".format(h1, mismatch)
                h2 = "{0} mismatch={1}".format(h2, mismatch)
            print("@{0}\n{1}\n+\n{2}".format(h1, s1, q1), file=outfiles[b.name][0])
            print("@{0}\n{1}\n+\n{2}".format(h2, s2, q2), file=outfiles[b.name][1])

    for key in outfiles:
        o1, o2 = outfiles[key]
        o1.close()
        o2.close()

    return total, counter

def demux(barcode_filename, r1, r2=None):
    barcode_collection = BarcodeCollection()
    barcode_collection.parse_tsv(barcode_filename)
    
    if r2 is None:
        total, counter = demux_single(barcode_collection, r1)
    else:
        total, counter = demux_pair(barcode_collection, r1, r2)

    expected = total / (len(counter)-1)
    thresh = 0.3
    above_expected = expected * (1 + thresh)
    below_expected = expected * (1 - thresh)
    print(total, "non-phiX reads processed")
    print("")
    print("sample", "count", "percent", "flag", sep='\t')
    for k, v in counter.most_common():
        tag = ""
        if v >= above_expected:
            tag = "above"
        elif v <= below_expected:
            tag = "below"
        print(k, v, v/total, tag, sep='\t')

    return total, counter
