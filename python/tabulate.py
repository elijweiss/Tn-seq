#!/usr/bin/env python
# Purpose: This program reads an input file listing positions
# and reads at each position for any number of Tn-Seq runs
#
# Copyright (c) 2014 University of Washington

#-----------------------------------------------------
# modules
#-----------------------------------------------------
import sys
import optparse
import os.path
import re
import collections
import common

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
EXCLUDE_3P = 5
EXCLUDE_5P = 90
OUT_SUFFIX = "_TabbyGn" + str(EXCLUDE_3P) + "-" + str(EXCLUDE_5P) + ".txt"
MIN_READS_PER_POS = 0
INCLUDE_YES = "1"

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-i", "--infile", action="store", type="string", dest="infile", help="path to annotated hits file")
    parser.add_option("-a", "--annofiles", action="store", type="string", dest="annofiles", help="path to .ptt annotation file(s) (comma-separated list if using more than one; order must match sequences in reference fasta)")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="path to output file with merged read counts (default <infile>" + OUT_SUFFIX + ")")
    parser.add_option("-f", "--fasta", action="store", type="string", dest="fasta", help="path to reference fasta containing replicon names (if using multiple replicons)")

    opts, args = parser.parse_args()
    if (opts.infile is None or opts.annofiles is None or opts.fasta is None):
        parser.print_help()
        exit(1)

    # Create the default output filename if not provided
    if (opts.outfile is None):
        parts = os.path.splitext(os.path.basename(opts.infile))
        opts.outfile = parts[0] + OUT_SUFFIX + parts[1]

    return opts, args

def read_hits_file(infile):
    hits = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(int)))
    allhits = collections.defaultdict(lambda: collections.defaultdict(int))
    with open(infile, "r") as fh:
        header = fh.readline()
        (replicon, eff_pos, direction, pid, locus, rel_pos, include, count_columns) = header.rstrip().split('\t', 7)
        run_names = count_columns.split('\t')
        num_runs = len(run_names)
        print str(num_runs) + " columns of read counts"

        count = 0
        for line in fh:
            (replicon, eff_pos, direction, pid, locus, rel_pos, include, count_columns) = line.rstrip().split('\t', 7)
            if include != INCLUDE_YES:
                continue
            readcounts = map(float, count_columns.split('\t'))

            # Only include insertions within genes - format #(#)
            mobj = re.match("^(\d+)\((\d+)\)", rel_pos)
            if mobj:
                # Acceptable window within ORF
                relpos_r = float(mobj.group(1)) / float(mobj.group(2))
                accept_pos = relpos_r * 100 >= EXCLUDE_3P and relpos_r * 100 <= EXCLUDE_5P
            else:
                continue

            count += 1
            hit_any_run = 0
            for col in range(0, num_runs):
                run = run_names[col]
                hit_any_run = 1
                if readcounts[col] <= 0 or readcounts[col] < MIN_READS_PER_POS:
                    continue
                # Hits anywhere within ORF
                hits[pid][run]['hits'] += 1
                hits[pid][run]['reads'] += readcounts[col]
                if accept_pos:
                    # Hits within acceptable window
                    hits[pid][run]['whits'] += 1
                    hits[pid][run]['wreads'] += readcounts[col]
            if hit_any_run:
                allhits[pid]['ahits'] += 1
                if accept_pos:
                    allhits[pid]['awhits'] += 1

    print "tabulated " + str(count) + " hits"
    return (run_names, hits, allhits)

def write_tabulated(outfile, runs, annotations, hits, allhits):
    count = 0
    with open(outfile, "w") as outfh:
        pct_range = str(EXCLUDE_3P) + "-" + str(EXCLUDE_5P)
        out_header = "Replicon\tLocus\tPID\tFrom\tTo\tStrand\tLength\tGene\tCode\tCog\tProduct\tHlocsORF_AllRuns\tHlocs" + pct_range + "_AllRuns"
        for run in runs:
            out_header += "\thTot_" + run + "\trTot_" + run + "\th" + pct_range + "_" + run + "\tr" + pct_range + "_" + run
        out_header += "\n"
        outfh.write(out_header)

        for replicon, annos in annotations.iteritems():
            for pid, pinfo in sorted(annos.iteritems(), key=lambda (k,v): int(v['startpos'])):
                outline = '\t'.join(map(str, [replicon, pinfo['locus_tag'], pid, pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['length'], pinfo['info'], allhits[pid]['ahits'], allhits[pid]['awhits']]))

                for run in runs:
                    outline += "\t" + str(hits[pid][run]['hits'])
                    outline += ("\t%.8f" % hits[pid][run]['reads']).rstrip('0').rstrip('.')
                    outline += "\t" + str(hits[pid][run]['whits'])
                    outline += ("\t%.8f" % hits[pid][run]['wreads']).rstrip('0').rstrip('.')

                outfh.write(outline + '\n')
                count += 1

    print "\rwrote " + str(count) + " positions to " + outfile
    return count

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()
    print "Tabulating " + opts.infile

    all_replicons = common.read_replicon_names(opts.fasta)

    # Separate multiple annotation files
    annofile_list = opts.annofiles.split(",")
    if len(annofile_list) != len(all_replicons):
        print "ERROR: wrong number of replicon names for replicon files"
        exit(1)
    print "Using annotation files " + str(annofile_list)

    annotations = common.read_annotations(annofile_list, all_replicons)

    (runs, hits, allhits) = read_hits_file(opts.infile)

    lines_written = write_tabulated(opts.outfile, runs, annotations, hits, allhits)

if __name__ == "__main__":
    main()
