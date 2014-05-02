#!/usr/bin/env python
# purpose: Read sum file and split columns into 2 additional files
#
# Copyright (c) 2014 University of Washington. All rights reserved.

#------------------------------------------------------------------------------
# system modules required
#------------------------------------------------------------------------------
import sys
import re

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def split_file(infile, outfile_reads, outfile_q0):
    lcount = 0
    ifh = open(infile, 'r')
    ofh_reads = open(outfile_reads, 'w')
    ofh_q0 = open(outfile_q0, 'w')
    for line in ifh:
        lcount += 1 
        if lcount == 1:
            outline_reads = "Replicon\tPosition of Insertion\tDirection\tReads(Norm)\n"
            outline_q0 = "Replicon\tPosition of Insertion\tDirection\tQ=0 Reads(Norm)\n"
        else:
            # Grap appropriate columns for output files
            (replicon, position, direction, reads, q0reads) = line.rstrip().split('\t')
            outline_reads = replicon + "\t" + position + "\t" + direction + "\t" + reads + "\n"
            outline_q0 = replicon + "\t" + position + "\t" + direction + "\t" + q0reads + "\n"

        ofh_reads.write(outline_reads)
        ofh_q0.write(outline_q0)

    print "Split " + str(lcount) + " lines"

#------------------------------------------------------------------------------
# main
#------------------------------------------------------------------------------
def main():
    if len(sys.argv) < 2 or len(sys.argv) == 3:
        print "Usage: " + sys.argv[0] + " <input file> [ <reads output file> <q=0 output file> ]"
        exit(1)

    infile = sys.argv[1]

    # Set output filenames to defaults if not provided
    p = re.compile("([\w\-\/]+)\.\w+")
    m = p.search(infile)
    prefix = m.group(1)
    if len(sys.argv) > 2:
        outfile_reads = sys.argv[2]
        outfile_q0 = sys.argv[3]
    else:
        outfile_reads = prefix + "_all.txt"
        outfile_q0 = prefix + "_q0.txt"

    print "Splitting reads summary into total reads and quality-0 reads"
    split_file(infile, outfile_reads, outfile_q0)

if __name__ == "__main__":
    main()
