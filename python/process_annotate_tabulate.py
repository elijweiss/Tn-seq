#!/usr/bin/env python
# Purpose: Run the annotation and tablulation steps. This can be run
# after several Tn-seq runs have been mapped using the mapping script.
# Takes any number of fastq files and produces 2 files of analysis results
#
# Copyright (c) 2014 University of Washington

#-----------------------------------------------------
# modules
#-----------------------------------------------------
import sys
import optparse
import os.path
import shutil
import multiprocessing
import common

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options] reads_list_1 reads_list_2 ..."
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-r", "--reffile", action="store", type="string", dest="reference_fa", help="path to reference genome fasta")
    parser.add_option("-a", "--annofiles", action="store", type="string", dest="annofiles", help="path to reference .ptt annotation file(s) (comma-separated list if using more than one; order must match sequences in reference fasta)")
    parser.add_option("-o", "--outfile_anno", action="store", type="string", dest="outfile_anno", default=common.OUTFILE_ANNO, help="path to final annotated output file (default: " + common.WORKING_DIR + "/" + common.OUTFILE_ANNO + ")")
    parser.add_option("-p", "--outfile_tab", action="store", type="string", dest="outfile_tab", default=common.OUTFILE_TAB, help="path to final counts tabulated by gene (default: " + common.WORKING_DIR + "/" + common.OUTFILE_TAB + ")")
    parser.add_option("-w", "--workingdir", action="store", default="", dest="workdir", help="working directory for input and output files (default: " + common.WORKING_DIR + ")")

    opts, args = parser.parse_args()
    if (opts.reference_fa is None or opts.annofiles is None):
        parser.print_help()
        exit(1)

    return opts, args

# Count the sequences in the fasta file
def count_fasta_seqs(fasta):
    count = 0
    with open(fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count

# Run the Tn-seq annotation pipeline
def process(infiles, opts):
    # Split summary files
    split_files_all = list()
    split_files_q0 = list()
    for norm_file in infiles:
        outfile_all = common.add_suffix(norm_file, common.ALL_SUFFIX)
        outfile_q0 = common.add_suffix(norm_file, common.Q0_SUFFIX)
        common.run_cmd(["python", "split_sum.py", norm_file, outfile_all, outfile_q0])
        split_files_all.append(outfile_all)
        split_files_q0.append(outfile_q0)

    # Compile summary sets
    readscomp = os.path.join(opts.workdir, common.READSCOMP)
    cmd = ["python", "compile_sets.py"]
    cmd.extend(split_files_all)
    cmd.extend(split_files_q0)
    cmd.append(readscomp)
    common.run_cmd(cmd)

    # Annotate positions
    reads_anno = common.add_suffix(readscomp, common.ANNO_SUFFIX)
    cmd = ["python", "annotate.py", "--infile", readscomp, "--annofiles", opts.annofiles, "--outfile", reads_anno, "--fasta", opts.reference_fa]
    common.run_cmd(cmd)

    # Merge annotations and counts
    anno_hits_file = os.path.join(opts.workdir, opts.outfile_anno)
    anno_reduced = common.add_suffix(reads_anno, common.TAB_SUFFIX)
    common.run_cmd(["python", "merge_anno.py", "--reads_file", readscomp, "--anno_file", reads_anno, "--outfile1", anno_hits_file, "--outfile2", anno_reduced])

    # Tabulate by gene
    tab_file = os.path.join(opts.workdir, opts.outfile_tab)
    cmd = ["python", "tabulate.py", "--infile", anno_reduced, "--annofiles", opts.annofiles, "--outfile", tab_file, "--fasta", opts.reference_fa]
    common.run_cmd(cmd)

    print "Processing complete. Output in " + anno_hits_file + " and " + tab_file


#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()
    infiles = args

    annofile_list = opts.annofiles.split(",")
    if count_fasta_seqs(opts.reference_fa) != len(annofile_list):
        print "Number of annotation files must match the number of sequences in the reference fasta"
        exit(1)

    process(infiles, opts)

if __name__ == "__main__":
    main()
