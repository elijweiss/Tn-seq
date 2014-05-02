#!/usr/bin/env python
# Purpose: Combine one or more sets of read counts into a single file.
# A set consists of two files, one with all reads and one with q=0 reads
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
SHOW_TOTALS = False

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def read_files(infiles):
    totals = dict()
    filereads = dict()
    for infile in infiles:
        print "reading file " + infile
        with open(infile, "r") as fh:
            try:
                header = fh.readline()
                for line in fh:
                    (replicon, pos, strand, readcount) = line.rstrip().split("\t")
                    totals[(replicon, pos, strand)] = totals.get((replicon, pos, strand), 0) + float(readcount)
                    filereads[(infile, replicon, pos, strand)] = readcount
            except StopIteration:
                break
    return (totals, filereads)

def write_compiled(totals, filereads, infiles, outfile):
    with open(outfile, "w") as fh:
        if SHOW_TOTALS:
            total_header = "\tAllReads"
        else:
            total_header = ""
        file_list = "\t".join(map(os.path.basename, infiles))
        fh.write("Replicon\tPosition\tDirection\t" + file_list + total_header + "\n")
        for (replicon, pos, strand) in sorted(totals, key=totals.get, reverse=True):
            outline = replicon + "\t" + pos + "\t" + strand
            for infile in infiles:
                numreads = filereads.get((infile, replicon, pos, strand), 0)
                outline += "\t" + str(numreads)
            if SHOW_TOTALS:
                total_val = "\t" + str(totals[(pos, strand)])
            else:
                total_val = ""
            fh.write(outline + total_val + "\n")
    print "total positions tabulated: " + str(len(totals))

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    if len(sys.argv) < 3:
        print "Usage: " + sys.argv[0] + " <input file(s)> <output file>"
        exit(1)
    infiles = sys.argv[1:-1]
    outfile = sys.argv[-1]

    print "Compiling sets of read counts"
    (totals, filereads) = read_files(infiles)

    write_compiled(totals, filereads, infiles, outfile)

if __name__ == "__main__":
    main()
