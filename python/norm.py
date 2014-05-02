#!/usr/bin/env python
# Purpose: Normalize read counts per location based on the
# total mapped reads for the set
#
# Copyright (c) 2014 University of Washington

#-----------------------------------------------------
# modules
#-----------------------------------------------------
import sys
import optparse
import os.path

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
OUT_SUFFIX = "_norm"
PROGRESS_FREQ = 10000
NORM_FACTOR = 10000000

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-i", "--infile", action="store", type="string", dest="infile", help="path to input file")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="path to normalized output file")
    parser.add_option("-n", "--norm-factor", action="store", type="int", dest="norm_factor", default=NORM_FACTOR, help="Normalization factor (default " + str(NORM_FACTOR) + ")")

    opts, args = parser.parse_args()
    if (opts.infile is None):
        parser.print_help()
        exit(1)

    # Create output filename if not provided
    if opts.outfile is None:
        parts = os.path.splitext(os.path.basename(opts.infile))
        opts.outfile = parts[0] + OUT_SUFFIX + parts[1]

    return opts, args

def normalize(infile, outfile, norm_factor):
    in_fh = open(infile, "r")
    out_fh = open(outfile, "w")
    count = 0
    total_reads = 0
    data = list()

    # First we need the total reads for this set
    for line in in_fh.readlines():
        if count == 0:
            header = line.rstrip()
        else:
            (replicon, position, direction, reads, q0) = line.rstrip().split("\t")
            total_reads += int(reads)
            data.append((replicon, position, direction, reads, q0))
        count += 1

    print "Total reads: " + str(total_reads)

    # Now normalize on the total
    count = 0
    out_fh.write(header + "\n")
    for (replicon, position, direction, reads, q0) in data:
        reads_n = float(reads) * int(norm_factor) / total_reads
        q0_n = float(q0) * int(norm_factor) / total_reads
        if q0_n == 0:
            q0_n = "0"
        normline = "\t".join((replicon, position, direction, str(reads_n), str(q0_n)))
        out_fh.write(normline + "\n")
        count += 1

    in_fh.close()
    out_fh.close()

    print "Normalized " + str(count) + " read counts"

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()

    print "Normalizing reads with normalization factor " + str(opts.norm_factor)
    normalize(opts.infile, opts.outfile,  opts.norm_factor)

if __name__ == "__main__":
    main()
