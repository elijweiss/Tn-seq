#!/usr/bin/env python
# Purpose: Filter a fastq file:
# Index tn seq must be no more than 1 base different from expected tn sequence.
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
FQ_LINES = 4
PROGRESS_FREQ = 10000
OUT_SUFFIX = "_iPass"
OTHER = "other"

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-e", "--end1", action="store", type="string", dest="e1file", help="path to first end fastq file")
    parser.add_option("-f", "--end2", action="store", type="string", dest="e2file", help="path to second end fastq file")
    parser.add_option("-i", "--index", action="store", type="string", dest="ifile", help="path to index fastq file")
    parser.add_option("-o", "--outfile1", action="store", type="string", dest="outfile1", help="path to output filtered 1st-end file (default <infile1>" + OUT_SUFFIX + ")")
    parser.add_option("-p", "--outfile2", action="store", type="string", dest="outfile2", help="path to output filtered 2nd-end file (default <infile2>" + OUT_SUFFIX + ")")
    parser.add_option("-s", "--tn_end", action="store", type="string", dest="tn_end", help="transposon end sequence")

    opts, args = parser.parse_args()
    if (opts.e1file is None or opts.ifile is None):
        parser.print_help()
        exit(1)

    # Create the default output filenames if not provided
    if (opts.outfile1 is None):
        parts = os.path.splitext(os.path.basename(opts.e1file))
        opts.outfile1 = parts[0] + OUT_SUFFIX + parts[1]
    if (opts.e2file and opts.outfile2 is None):
        parts = os.path.splitext(os.path.basename(opts.e2file))
        opts.outfile2 = parts[0] + OUT_SUFFIX + parts[1]

    return opts, args

def filter_tnend(e1file, e2file, ifile, outfile1, outfile2, tnend):
    tn_end_odd_length = 0
    ind_end_odd_length = 0
    similar_seqs = get_similar_seqs(tnend)
    similar_counts = dict()
    e1_match_counts = dict()
    e1_match_total = 0
    readcount = 0
    i_ok_count = 0
    index_length = 0
    end_length = len(tnend)

    e1fh = open(e1file, "r")
    ifh = open(ifile, "r")
    o1fh = open(outfile1, "w")
    if e2file:
        e2fh = open(e2file, "r")
        o2fh = open(outfile2, "w")
    while(True):
        try:
            e1read = [e1fh.next() for x in range(FQ_LINES)]
            iread = [ifh.next() for x in range(FQ_LINES)]
            if e2file:
                e2read = [e2fh.next() for x in range(FQ_LINES)]
            if readcount == 0:
                index_length = len(iread[1].rstrip())
            read_tn_end = iread[1].rstrip()[:end_length].upper()
            rest = iread[1].rstrip()[end_length:].upper()

            if len(read_tn_end) != end_length:
                print "unusual tn end length: " + read_tn_end
                tn_end_odd_length += 1
            if len(rest) != index_length - end_length:
                print "unusual index read end length: " + rest
                ind_end_odd_length += 1

            if read_tn_end in similar_seqs:
                similar_counts[read_tn_end] = similar_counts.get(read_tn_end, 0) + 1
                i_ok_count += 1
                if re.match(rest, e1read[1].upper()):
                    e1_match_counts[read_tn_end] = e1_match_counts.get(read_tn_end, 0) + 1
                    e1_match_total += 1
                    o1fh.writelines(e1read)
                    if e2file:
                        o2fh.writelines(e2read)
            else:
                similar_counts['other'] = similar_counts.get('other', 0) + 1
                if re.match(rest, e1read[1].upper()):
                    e1_match_counts['other'] = e1_match_counts.get('other', 0) + 1

            readcount += 1
            if readcount % PROGRESS_FREQ == 0:
                sys.stdout.write("processed " + str(readcount) + " reads\r")
                sys.stdout.flush()

        except StopIteration:
            break

    e1fh.close()
    ifh.close()
    o1fh.close()
    if e2file:
        e2fh.close()
        o2fh.close()

    # Print report
    print "Total reads processed:  " + str(readcount)
    print "Total with 1st %d of i read matching acceptable sequence(s):  %d  (%.3f%%)" % (end_length, i_ok_count, (100.0 * i_ok_count / readcount))
    print "Total passing filter (e1-i match):  %d  (%.3f%%)" % (e1_match_total, (100.0 * e1_match_total / i_ok_count))
    print " Reads parsed for tn end sequences:"
    for seq in similar_seqs:
        if seq != OTHER:
            print_report_line(seq, similar_counts.get(seq, 0), e1_match_counts.get(seq, 0))
    print_report_line(OTHER, similar_counts.get(OTHER, 0), e1_match_counts.get(OTHER, 0))
    print "\nNumber with unusual length for end of tn read: " + str(tn_end_odd_length)
    print "Number with unusual length for end of index read: " + str(ind_end_odd_length)
    print "Wrote " + str(e1_match_total) + " filtered reads"

def print_report_line(seq, total_reads, e1_match_counts):
    if (total_reads) > 0:
        pct = "(%3.3f%%)" % (100.0 * e1_match_counts / total_reads)
    else:
        pct = ""
    print "  %s\t%12d\t%12d  %s" % (seq, total_reads, e1_match_counts, pct)


# Return a list of all strings that differ from the input string by at most one base
def get_similar_seqs(tnend):
    seqs = [ tnend ]
    seqs.extend([ tnend[:p] + b + tnend[p+1:] for p in range(len(tnend)) for b in "ACGTN" if tnend[p] != b ])
    return seqs

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()

    print "Filtering for Tn end sequence"
    filter_tnend(opts.e1file, opts.e2file, opts.ifile, opts.outfile1, opts.outfile2, opts.tn_end)

if __name__ == "__main__":
    main()
