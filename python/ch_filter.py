#!/usr/bin/env python
# Purpose: Filter a fastq file for chastity
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
OUT_SUFFIX = "ch"
FQ_LINES = 4
PROGRESS_FREQ = 10000

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-e", "--end1", action="store", type="string", dest="e1_file", help="path to first end fastq file")
    parser.add_option("-f", "--end2", action="store", type="string", dest="e2_file", help="path to second end fastq file")
    parser.add_option("-i", "--index", action="store", type="string", dest="index_file", help="path to index fastq file")
    parser.add_option("-o", "--outfile_e1", action="store", type="string", dest="outfile_e1", help="path to output filtered 1st-end file")
    parser.add_option("-p", "--outfile_e2", action="store", type="string", dest="outfile_e2", help="path to output filtered 2nd-end file")
    parser.add_option("-q", "--outfile_i", action="store", type="string", dest="outfile_i", help="path to output filtered index file")

    opts, args = parser.parse_args()
    if (opts.e1_file is None):
        parser.print_help()
        exit(1)

    # Create output filenames if not provided
    if (opts.outfile_e1 is None):
        parts = os.path.splitext(opts.e1_file)
        opts.outfile_e1 = parts[0] + OUT_SUFFIX + parts[1]
    if (opts.outfile_e2 is None and opts.e2_file is not None):
        parts = os.path.splitext(opts.e2_file)
        opts.outfile_e2 = parts[0] + OUT_SUFFIX + parts[1]
    if (opts.outfile_i is None and opts.index_file is not None):
        parts = os.path.splitext(opts.index_file)
        opts.outfile_i = parts[0] + OUT_SUFFIX + parts[1]

    return opts, args

def filter_chastity(e1_file, e2_file, index_file, e1_chaste, e2_chaste, index_chaste):
    e1_f = open(e1_file, "r")
    e1_chf = open(e1_chaste, "w")
    if (e2_file is not None):
        e2_f = open(e2_file, "r")
        e2_chf = open(e2_chaste, "w")
    if (index_file is not None):
        index_f = open(index_file, "r")
        index_chf = open(index_chaste, "w")

    reads = 0
    ch_all = 0
    ch_some = 0
    mismatch_names = 0

    while True:
        e1_lines = list()
        e2_lines = list()
        index_lines = list()
        for _ in range(0, FQ_LINES):
            e1_lines.append(e1_f.readline())
            if (e2_file is not None):
                e2_lines.append(e2_f.readline())
            if (index_file is not None):
                index_lines.append(index_f.readline())

        if len(e1_lines[0]) == 0:
            break

        (e1_name, ch_string) = e1_lines[0].rsplit(":", 1)
        e1_chaste = not ch_string.startswith("0")
        if (e2_file is None):
            e2_chaste = False
        else:
            (e2_name, ch_string) = e2_lines[0].rsplit(":", 1)
            e2_chaste = not ch_string.startswith("0")
        if (index_file is None):
            i_chaste = False
        else:
            (i_name, ch_string) = index_lines[0].rsplit(":", 1)
            i_chaste = not ch_string.startswith("0")

        if ((e2_file is not None and e1_name != e2_name) or
            (index_file is not None and e1_name != i_name)):
            mismatch_names += 1
            next()

        if e1_chaste and e2_chaste and i_chaste:
            ch_all += 1
        elif e1_chaste or e2_chaste or i_chaste:
            ch_some += 1

        if e1_chaste or e2_chaste or i_chaste:
            e1_chf.writelines(e1_lines)
            if (e2_file is not None):
                e2_chf.writelines(e2_lines)
            if (index_file is not None):
                index_chf.writelines(index_lines)

        reads += 1
        if reads % PROGRESS_FREQ == 0:
            sys.stdout.write("\r" + str(reads) + " filtered (" + str(ch_all + ch_some) + " chaste)")
            sys.stdout.flush()

    e1_f.close()
    e1_chf.close()
    if (e2_file is not None):
        e2_f.close()
        e2_chf.close()
    if (index_file is not None):
        index_f.close()
        index_chf.close()

    print "\nTotal reads analyzed: " + str(reads)
    print " Chaste reads: " + str(ch_all + ch_some)
    print "    Number chaste in all reads: " + str(ch_all)
    print "    Number chaste in only one: " + str(ch_some)
    print " Reads with mismatched names (discarded): " + str(mismatch_names)

    return ch_all + ch_some

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()

    if (opts.index_file is None):
        print "Filtering %s for chastity" % (opts.e1_file)
    else:
        print "Filtering %s and %s for chastity" % (opts.e1_file, opts.index_file)

    chaste_reads = filter_chastity(opts.e1_file, opts.e2_file, opts.index_file, opts.outfile_e1, opts.outfile_e2, opts.outfile_i)
    if chaste_reads == 0:
        print "ERROR: no chaste reads found"
        exit(1)

if __name__ == "__main__":
    main()
