#!/usr/bin/env python
# Purpose: Filter a fastq file:
# Read1 tn seq must be no more than 1 base different from expected tn sequence.
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
FQ_LINES = 4
PROGRESS_FREQ = 10000
OUT_SUFFIX = "_trim"
TN_END_SEQ = "AGACAG"
HAMM = 1

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-e", "--end1", action="store", type="string", dest="e1file", help="path to first end fastq file")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="path to output filtered fastq file (default <infile>" + OUT_SUFFIX + ")")
    parser.add_option("-s", "--seq", action="store", type="string", dest="tn_end_seq", help="transposon end sequence (default " + TN_END_SEQ + ")")

    opts, args = parser.parse_args()
    if (opts.e1file is None):
        parser.print_help()
        exit(1)

    if (opts.tn_end_seq is None):
        opts.tn_end_seq = TN_END_SEQ

    # Create the default output filename if not provided
    if (opts.outfile is None):
        parts = os.path.splitext(os.path.basename(opts.e1file))
        opts.outfile = parts[0] + OUT_SUFFIX + parts[1]

    return opts, args

def hamm_dist(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Sequences are of unequal length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def filter_reads(e1file, outfile, tn_end_seq):
    readcount = 0
    okcount = 0
    end_length = len(tn_end_seq)

    e1fh = open(e1file, "r")
    ofh = open(outfile, "w")
    while(True):
        try:
            # Filter out reads with too many differences. Trim end off before writing.
            e1read = [e1fh.next() for x in range(FQ_LINES)]
            readcount += 1
            read_end = e1read[1][:end_length].upper()
            if hamm_dist(read_end, tn_end_seq) <= HAMM:
                okcount += 1
                e1read[1] = e1read[1][end_length:]
                e1read[3] = e1read[3][end_length:]
                for line in e1read:
                    ofh.write(line)

            if readcount % PROGRESS_FREQ == 0:
                sys.stdout.write("\r" + str(readcount) + " processed")
                sys.stdout.flush()

        except StopIteration:
            break

    e1fh.close()
    ofh.close()

    # Print report
    print "\rTotal reads processed:  " + str(readcount)
    print "Total with 1st %d of read matching tn end seq:  %d  (%.3f%%)" % (end_length, okcount, (100.0 * okcount / readcount))
    print "Transposon ends trimmed"

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()

    print "Filtering for Tn end sequence. Expected sequence: " + opts.tn_end_seq
    filter_reads(opts.e1file, opts.outfile, opts.tn_end_seq)

if __name__ == "__main__":
    main()
