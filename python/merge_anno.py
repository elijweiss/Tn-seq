#!/usr/bin/env python
# Purpose: Merge read counts file with annotated file
# (assumes both files are in the same order with one line per insertion site)
#
# Copyright (c) 2014 University of Washington. All rights reserved.

#------------------------------------------------------------------------------
# system modules required
#------------------------------------------------------------------------------
import sys
import optparse
import os

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
OUT1_SUFFIX = "_rcomp_merged.csv"
OUT2_SUFFIX = "_HitsToTab.txt"
INCLUDE = 1

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-r", "--reads_file", action="store", type="string", dest="reads_file", help="Path to input file containing read counts")
    parser.add_option("-a", "--anno_file", action="store", type="string", dest="anno_file", help="Path to annotated read count file")
    parser.add_option("-o", "--outfile1", action="store", type="string", dest="outfile1", help="Merged output file")
    parser.add_option("-p", "--outfile2", action="store", type="string", dest="outfile2", help="Output for Tabulate script")

    opts, args = parser.parse_args()
    if (opts.reads_file is None or opts.anno_file is None):
        parser.print_help()
        exit(1)

    return opts, args

def merge_files(reads_file, anno_file, outfile1, outfile2):
    try:
        # Open input files and skip header lines
        reads_fh = open(reads_file, "r")
        reads_header = reads_fh.readline()
        anno_fh = open(anno_file, "r")
        firstline = anno_fh.readline()

        (replicon, position, direction, reads_columns) = reads_header.rstrip().split('\t', 3)
        reads_headings = '\t'.join(reads_columns.split('\t'))

        data = []
        ambig_data = []
        mismatches = 0
        for line_reads in reads_fh:
            line_anno = anno_fh.readline()
            (replicon, position, direction, reads_columns) = line_reads.rstrip('\n').split('\t', 3)
            (repl_anno, pos_anno, dir_anno, eff_pos, pid, locus, from_pos, to_pos, strand, gene, code, cog, product, rel_pos, notes) = line_anno.rstrip('\n').split('\t')

            # Save hits within multiple possible genes
            if eff_pos == "":
                outline1 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (position, eff_pos, direction, pid, locus, strand, gene, product, rel_pos, notes, reads_columns)
                ambig_data.append((replicon, position, direction, outline1))
                continue

            # Make sure the input files match at this line
            if replicon != repl_anno or position != pos_anno or direction != dir_anno:
                mismatches += 1
                continue

            outline1 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (position, eff_pos, direction, pid, locus, strand, gene, product, rel_pos, notes, reads_columns)
            outline2 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (replicon, eff_pos, direction, pid, locus, rel_pos, INCLUDE, reads_columns)
            data.append((replicon, eff_pos, direction, outline1, outline2))

        print "read " + str(len(data)) + " lines from each input file"
        print str(mismatches) + " mismatched lines found"

        out1_fh = open(outfile1, "w")
        out2_fh = open(outfile2, "w")
        out1_fh.write("Position of insertion(t+1)\tEffective Position\tDirection\tPID\tLocus\tStrand\tGene\tProduct\tPos Relative to Locus\tNotes\t" + reads_headings + "\n")
        out2_fh.write("Replicon\tEffective Position\tDirection\tPID\tLocus\tPosition relative to locus\tInclude\t" + reads_headings + "\n")

        data_sorted = sorted(data, key = lambda x: (x[0], int(x[1])))
        for (repl, eff_pos, direction, outline1, outline2) in data_sorted:
            out1_fh.write(outline1)
            out2_fh.write(outline2)

        # Write hits with multiple possible genes at the end
        ambig_data_sorted = sorted(ambig_data, key = lambda x: (x[0], int(x[1])))
        for (repl, position, direction, outline1) in ambig_data_sorted:
            out1_fh.write(outline1)

        print "wrote " + outfile1 + " and " + outfile2

    except IOError as e:
        print e

#------------------------------------------------------------------------------
# main
#------------------------------------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()

    # Create output filenames if not supplied
    if (opts.outfile1 is None):
        parts = os.path.splitext(os.path.basename(opts.reads_file))
        opts.outfile1 = parts[0] + OUT1_SUFFIX
    if (opts.outfile2 is None):
        parts = os.path.splitext(os.path.basename(opts.reads_file))
        opts.outfile2 = parts[0] + OUT2_SUFFIX

    print "Merging read counts with annotations"
    merge_files(opts.reads_file, opts.anno_file, opts.outfile1, opts.outfile2)

if __name__ == "__main__":
    main()
