#!/usr/bin/env python
# Purpose: Reads mapping file (.sam format) and produces a file containing
# read counts per insertion location
#
# Copyright (c) 2014 University of Washington

#-----------------------------------------------------
# modules
#-----------------------------------------------------
import sys
import optparse
import os.path
import operator

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
FQ_LINES = 4
PROGRESS_FREQ = 10000
OUT_SUFFIX = "_sum.txt"
SWAP_DIR = False

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-i", "--infile", action="store", type="string", dest="infile", help="path to input sam file")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="path to output counts file (default <infile>" + OUT_SUFFIX + ")")

    opts, args = parser.parse_args()
    if (opts.infile is None):
        parser.print_help()
        exit(1)

    # Create the default output filename if not provided
    if (opts.outfile is None):
        parts = os.path.splitext(os.path.basename(opts.infile))
        opts.outfile = parts[0] + OUT_SUFFIX

    return opts, args

# Parse .sam file and count mapped reads per position
def count_reads(infile):
    with open(infile, "r") as fh:
        count = 0
        readlength = 0
        poscounts = dict()
        zerocounts = dict()
        for line in fh:
            if line.startswith("@"):
                continue
            (qname, flag, rname, pos, mapq, cigar, mrnm, mpos, isize, seq, qual) = line.rstrip().split("\t", 10)
            # Sequencing of Tn is usually 5' to 3', so a + mapping is considered Reverse
            if flag == "4":
                continue
            if flag == "16":
                strand = "F"
            else:
                strand = "R"
            replicon = rname
            if not readlength:
                readlength = len(seq)
            # Adjust to the Tn start position for F reads
            if strand == "F":
                pos = str(int(pos) + (readlength - 1))
            if SWAP_DIR:
                if strand == "F":
                    strand = "R"
                else:
                    strand = "F"
            poscounts[(replicon, pos, strand)] = poscounts.get((replicon, pos, strand), 0) + 1
            if int(mapq) == 0:
                zerocounts[(replicon, pos, strand)] = zerocounts.get((replicon, pos, strand), 0) + 1
            count += 1
            if count % PROGRESS_FREQ == 0:
                sys.stdout.write("read " + str(count) + " mappings\r")
                sys.stdout.flush()
    print "Total mapped reads analyzed: " + str(count)
    return (poscounts, zerocounts)

def print_counts(outfile, poscounts, zerocounts):
    with open(outfile, "w") as out:
        out.write("Replicon\tPosition\tDirection\tTotal Reads\tQ=0 Reads\n")
        count = 0
        sortedkeys = sorted(poscounts.keys(), key=lambda x: x[2], reverse=True)
        sortedkeys.sort(key=lambda x: int(x[1]))
        sortedkeys.sort(key=lambda x: x[0])
        for (replicon, pos, strand) in sortedkeys:
            totalreads = poscounts[(replicon, pos, strand)]
            q0reads = zerocounts.get((replicon, pos, strand), 0)
            out.write(replicon + "\t" + pos + "\t" + strand + "\t" + str(totalreads) + "\t" + str(q0reads) + "\n")
            count += 1
        print "Positions written: " + str(count)

def revcomp(seq):
    complements = string.maketrans("acgtACGT", "tgcaTGCA")
    rc = seq.translate(complements)[::-1]
    return rc

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()

    print "Summarizing mappings"
    (poscounts, zerocounts) = count_reads(opts.infile)
    print_counts(opts.outfile, poscounts, zerocounts)

    if len(poscounts) == 0:
        print "ERROR: no mapped reads found in input file"
        exit(1)

if __name__ == "__main__":
    main()
