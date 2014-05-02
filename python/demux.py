#!/usr/bin/env python
# Purpose: De-multiplex reads into multiple files
# Output filenames will contain demultiplexing barcode:
# if input file is reads.fq, output will be in the form
# reads_AGACAG.fq.
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
    parser.add_option("-s", "--seqs", action="store", type="string", dest="seqs_file", help="path to list of sample index sequences")
    parser.add_option("-o", "--out_prefix", action="store", type="string", dest="out_prefix", help="optional prefix for output files")

    opts, args = parser.parse_args()
    if (opts.e1_file is None or (opts.e2_file is None and opts.index_file is None)):
        parser.print_help()
        exit(1)

    return opts, args

# Read in the list of possible index sequences for de-multiplexing
def read_seqs(seqs_file):
    seqlist = list()
    with open(seqs_file, "r") as fh:
        for line in fh:
            seq = line.rstrip()
            seqlist.append(seq)
    return seqlist

def hamm_dist(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Sequences are of unequal length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# Find the barcode in the list that matches the input sequence
def find_best_seq(sequence, seqlist):
    best_seq = None
    for candidate in seqlist:
        if (hamm_dist(candidate, sequence) <= 1):
            best_seq = candidate
            break
    return best_seq

# Split the input into multiple output files based on the index sequences
def demultiplex(e1_file, e2_file, index_file, seqlist):
    # Open input files appropriate to script mode
    e1_fh = open(e1_file, "r")
    keep_index = False
    if (e2_file is None):
        # Demux using index read, output 1st-end reads
        print "De-multiplexing using index reads"
        d_fh = open(index_file, "r")
    elif (index_file is None):
        # Demux using 2nd-end read, output 1st-end reads
        print "De-multiplexing using 2nd-end reads"
        d_fh = open (e2_file, "r")
    elif (index_file is not None and e2_file is not None):
        # Demux using 2nd-end read, output 1st-end and index read
        print "De-multiplexing using 2nd-end reads"
        d_fh = open (e2_file, "r")
        i_fh = open(index_file, "r")
        keep_index = True

    # Open output files with appropriate names
    e1_outfiles = dict()
    i_outfiles = dict()
    (out_prefix_e1, out_extn) = os.path.splitext(e1_file)
    if index_file is not None:
        (out_prefix_i, _) = os.path.splitext(index_file)
    for seq in seqlist:
        outfile_e1 = out_prefix_e1 + "_" + seq + out_extn
        e1_outfiles[seq] = open(outfile_e1, "w")
        if (index_file is not None and e2_file is not None):
            outfile_i = out_prefix_i + "_" + seq + out_extn
            i_outfiles[seq] = open(outfile_i, "w")

    counts = dict()
    count_total = 0
    count_rejects = 0
    barcode_length = len(seqlist[0])
    while True:
        d_lines = list()
        e1_lines = list()
        i_lines = list()
        for _ in range(0, FQ_LINES):
            e1_lines.append(e1_fh.readline())
            d_lines.append(d_fh.readline())
            if e2_file is None:
                i_lines.append(d_lines[-1])
            else:
                if (keep_index):
                    i_lines.append(i_fh.readline())
        if len(d_lines[0]) == 0:
            break
        count_total += 1

        read_seq = d_lines[1].rstrip()
        bin_seq = find_best_seq(read_seq[:barcode_length], seqlist)
        if bin_seq is None:
            count_rejects += 1
            continue

        e1_outfiles[bin_seq].writelines(e1_lines)
        if keep_index:
            i_outfiles[bin_seq].writelines(i_lines)

        counts[bin_seq] = counts.get(bin_seq, 0) + 1
        if count_total % PROGRESS_FREQ == 0:
            sys.stdout.write("\r" + str(count_total) + " reads processed")
            sys.stdout.flush()

    e1_fh.close()
    d_fh.close()
    for seq in seqlist:
        e1_outfiles[seq].close()
        if (index_file is not None and e2_file is not None):
            i_outfiles[seq].close()
    if keep_index:
        i_fh.close()

    print "\nTotal reads processed: " + str(count_total)
    print "Reads with invalid index reads: " + str(count_rejects)
    print "Counts by file:"
    for seq in seqlist:
        print "  " + e1_outfiles[seq].name + ": " + str(counts[seq])

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()
    
    seqlist = read_seqs(opts.seqs_file)
    demultiplex(opts.e1_file, opts.e2_file, opts.index_file, seqlist)

if __name__ == "__main__":
    main()
