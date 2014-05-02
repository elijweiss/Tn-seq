#!/usr/bin/env python
# Purpose: Count occurrences of the first n bases of the index reads.
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
FQ_LINES = 4
PROGRESS_FREQ = 100000
TN_END_LENGTH = 6
OUT_SUFFIX = "_index_log.xls"

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-i", "--infile", action="store", type="string", dest="infile", help="path to input fastq file")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="path to output counts file (default <infile>" + OUT_SUFFIX + ")")
    parser.add_option("-l", "--tn_end_length", action="store", type="int", default=TN_END_LENGTH, dest="tn_end_length", help="number of bases in index (default " + str(TN_END_LENGTH) + ")")

    opts, args = parser.parse_args()
    if (opts.infile is None):
        parser.print_help()
        exit(1)

    # Create the default output filename if not provided
    if (opts.outfile is None):
        parts = os.path.splitext(os.path.basename(opts.infile))
        opts.outfile = parts[0] + OUT_SUFFIX

    return opts, args

# Look at the first n bases of each read in the fastq and keep a tally
# of occurrences of each sequence
def count_indexes(infile, tn_end_length):
    with open(infile, "r") as fh:
        bins = dict()
        total = 0
        while(True):
            try:
                read = [fh.next() for x in xrange(FQ_LINES)]
                first_bases = read[1][:tn_end_length]
                bins[first_bases] = bins.get(first_bases, 0) + 1
                total += 1
                if total % PROGRESS_FREQ == 0:
                    sys.stdout.write("\rCounted " + str(total) + " reads")
                    sys.stdout.flush()
            except StopIteration:
                break
        print "\rCounted " + str(total) + " reads"
        return bins

def write_counts(bins, outfile):
    with open(outfile, "w") as out:
        for index in sorted(bins, key=bins.get, reverse=True):
            out.write(index + "\t" + str(bins[index]) + "\n")
        print "Wrote " + str(len(bins)) + " counts to " + outfile

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()

    print "Counting sequences of index reads"
    bins = count_indexes(opts.infile, opts.tn_end_length)
    write_counts(bins, opts.outfile)

if __name__ == "__main__":
    main()
