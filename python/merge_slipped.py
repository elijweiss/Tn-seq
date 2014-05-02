#!/usr/bin/env python
# Purpose: Merge "slipped" reads (adjacent, cooriented locations with
# highly disproportionate read counts)
# The rest of the index must match the first bases of the first-end read.
#
# Copyright (c) 2014 University of Washington

#-----------------------------------------------------
# modules
#-----------------------------------------------------
import sys
import optparse
import os.path
import re

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
PROGRESS_FREQ = 10000
OUT_SUFFIX = "_merge.txt"
PROXIMITY = 2
PROPORTION = 0.04

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-i", "--infile", action="store", type="string", dest="infile", help="path to read counts summary file")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="path to output file with merged read counts (default <infile>" + OUT_SUFFIX + ")")

    opts, args = parser.parse_args()
    if (opts.infile is None):
        parser.print_help()
        exit(1)

    # Create the default output filename if not provided
    if (opts.outfile is None):
        parts = os.path.splitext(os.path.basename(opts.infile))
        opts.outfile = parts[0] + OUT_SUFFIX + parts[1]

    return opts, args

def read_summary(infile):
    poscounts = dict()
    with open(infile, "r") as fh:
        try:
            header = fh.readline()
            for line in fh:
                (replicon, pos, strand, reads, q0reads) = line.rstrip().split('\t')
                poscounts[(replicon, int(pos), strand)] = (int(reads), int(q0reads))
        except StopIteration:
            pass
    return (poscounts, header)

def merge_slips(poscounts):
    prevrepl = None
    prevpos = None
    prevstrand = None
    prevreads = 0
    prevq0 = 0
    for (replicon, pos, strand) in sorted(poscounts, key=lambda x: int(x[1])):
        (reads, q0reads) = poscounts[(replicon, pos, strand)]
        if prevrepl is not None:
            (prevreads, prevq0) = poscounts[(prevrepl, prevpos, prevstrand)]

        # We should merge if the groups of reads are close enough together,
        # on the same strand, and disproportionate enough.
        # Merge the smaller count into the larger.
        if (replicon == prevrepl and strand == prevstrand and abs(pos - prevpos) <= PROXIMITY
            and ( float(reads) / prevreads <= PROPORTION or float(prevreads) / reads <= PROPORTION )):
            if reads > prevreads:
                poscounts[(replicon, pos, strand)] = ( (reads + prevreads), (q0reads + prevq0) )
                poscounts.pop( (prevrepl, prevpos, prevstrand) )
                prevrepl = replicon
                prevpos = pos
                prevstrand = strand
            else:
                poscounts[(prevrepl, prevpos, prevstrand)] = ( (reads + prevreads), (q0reads + prevq0) )
                poscounts.pop( (replicon, pos, strand) )
        else:
            prevrepl = replicon
            prevpos = pos
            prevstrand = strand
    return poscounts

def write_merged(mgcounts, header, outfile):
    with open(outfile, "w") as fh:
        fh.write(header)
        for (replicon, pos, strand) in sorted(mgcounts, key=lambda x: int(x[1])):
            (reads, q0reads) = mgcounts[(replicon, pos, strand)]
            fh.write('\t'.join( (replicon, str(pos), strand, str(reads), str(q0reads)) ) + '\n')

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()
    print "Merging slipped reads\nProximity setting: " + str(PROXIMITY) + "\nProportion setting: " + str(PROPORTION)

    (poscounts, header) = read_summary(opts.infile)
    print str(len(poscounts)) + " positions in input file"
    mgcounts = merge_slips(poscounts)
    write_merged(mgcounts, header, opts.outfile)
    print str(len(mgcounts)) + " positions written to " + opts.outfile

if __name__ == "__main__":
    main()
