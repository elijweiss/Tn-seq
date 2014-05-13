#!/usr/bin/env python
# Purpose: Run the Tn-seq mapping pipeline. Produces files
# listing reads per location. The output files can be fed
# to the annotation script for analysis.
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
    usage = "usage: %prog [options] firstend_fastq_1 index_fastq_1 secondend_fastq_1 ..."
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-r", "--reffile", action="store", type="string", dest="reference_fa", help="path to reference genome fasta")
    parser.add_option("-j", "--tn_verify_by_read1", action="store_true", default=False, dest="verify_r1", help="use 1st-end read to verify transposon end (default: False)")
    parser.add_option("-i", "--tn_verify_by_index", action="store_true", default=False, dest="verify_i", help="use index read to verify transposon end (default: False)")
    parser.add_option("-t", "--tn_end", action="store", type="string", dest="tn_end_seq", help="expected transposon end sequence")
    parser.add_option("-d", "--demux_index", action="store_true", default=False, dest="demux_i", help="use index read to demultiplex (default: False)")
    parser.add_option("-e", "--demux_read2", action="store_true", default=False, dest="demux_r2", help="use 2nd end read to demultiplex (default: False)")
    parser.add_option("-b", "--barcodefile", action="store", type="string", dest="barcodefile", help="path to file listing expected barcode sequences")
    parser.add_option("-c", "--chastity", action="store_true", default=False, dest="dochastity", help="run chastity filter (default: False)")
    parser.add_option("-n", "--normfactor", action="store", type="int", default=common.NORM_FACTOR, dest="norm_factor", help="read count normalization factor (default: " + str(common.NORM_FACTOR) + ") (0 = don't normalize)")
    parser.add_option("-s", "--merge_slipped", action="store_true", default=False, dest="merge_slipped", help="merge slipped reads (default: False)")
    parser.add_option("-u", "--use_bowtie", action="store_true", default=False, dest="use_bowtie", help="map reads using Bowtie (default: use BWA)")
    parser.add_option("-w", "--workingdir", action="store", default=common.WORKING_DIR, dest="workdir", help="working directory for input and output files (default: " + common.WORKING_DIR + ")")

    opts, args = parser.parse_args()
    if (opts.reference_fa is None):
        parser.print_help()
        exit(1)

    return opts, args

def check_hash(hash_log, expected_seq):
    # Check that most of the index sequences match the expected sequence
    with open(hash_log) as fh:
        count_expected = 0
        count_rest = 0
        for line in fh:
            (index_seq, count) = line.rstrip().split('\t')
            if index_seq == expected_seq:
                count_expected = count
            else:
                count_rest += int(count)
        index_seq_ratio = float(count_expected) / (float(count_expected) + float(count_rest))
        return ((index_seq_ratio >= common.IDX_RATIO), index_seq_ratio)

# Read in the list of possible index sequences for de-multiplexing
def read_barcodes(barcodefile):
    seqlist = list()
    with open(barcodefile, "r") as fh:
        for line in fh:
            seq = line.rstrip()
            seqlist.append(seq)
    return seqlist

# Check fasta headers: BWA parses headers containing pipes and spaces so we recommend using simpler headers
def headers_ok(fasta):
    ok = True
    replicons = common.read_replicon_names(fasta)
    for replicon in replicons.values():
        if "|" in replicon:
            ok = False
    return ok

# Run the Tn-seq mapping pipeline
def process(infiles, opts):
    scriptpath = os.path.dirname(os.path.abspath(__file__))
    label = os.path.splitext(os.path.basename(infiles[0]))[0]

    # Group the input files (read1, index, read2) if using index reads and/or 2nd-end reads
    infile_groups = list()
    try:
        while len(infiles) > 0:
            fq1 = infiles.pop(0)
            group = (fq1, None, None)
            if opts.demux_r2 and (opts.demux_i or opts.verify_i):
                group = (fq1, infiles.pop(0), infiles.pop(0))
            elif opts.demux_r2:
                group = (fq1, None, infiles.pop(0))
            elif opts.demux_i or opts.verify_i:
                group = (fq1, infiles.pop(0), None)
            infile_groups.append(group)
    except IndexError:
        print "ERROR: wrong number of input files"
    print "Processing " + str(len(infile_groups)) + " Tn-seq fileset(s)"

    # Validation: make sure we have the required input files for the chosen options
    (r1_file, i_file, r2_file) = infile_groups[0]
    if (opts.verify_i or opts.demux_i) and i_file is None:
        print "Please supply index fastq files when using demultiplex-by-index or verify-by-index options"
        exit(1)
    if (opts.demux_r2) and r2_file is None:
        print "Please supply 2nd-end fastq files when using demultiplex-by-read2 option"
        exit(1)

    # Create the working directory if needed
    try:
        os.makedirs(opts.workdir)
    except OSError:
        if not os.path.isdir(opts.workdir):
            raise

    # Chastity filter
    chaste_files = list()
    if opts.dochastity:
        for (r1_file, i_file, r2_file) in infile_groups:
            fname = common.add_suffix(os.path.basename(r1_file), common.CHASTE_SUFFIX)
            r1_out = os.path.join(opts.workdir, fname)
            cmd = ["python", os.path.join(scriptpath, "ch_filter.py"), "--end1", r1_file, "--outfile_e1", r1_out]
            if r2_file is None:
                r2_out = None
            else:
                fname = common.add_suffix(os.path.basename(r2_file), common.CHASTE_SUFFIX)
                r2_out = os.path.join(opts.workdir, fname)
                cmd.extend(["--end2", r2_file, "--outfile_e2", r2_out])
            if i_file is None:
                i_out = None
            else:
                fname = common.add_suffix(os.path.basename(i_file), common.CHASTE_SUFFIX)
                i_out = os.path.join(opts.workdir, fname)
                cmd.extend(["--index", i_file, "--outfile_i", i_out])
            common.run_cmd(cmd)
            chaste_group = (r1_out, i_out, r2_out)
            chaste_files.append(chaste_group)
    else:
        for (r1_file, i_file, r2_file) in infile_groups:
            if r1_file:
                r1_out = os.path.join(opts.workdir, os.path.basename(r1_file))
                shutil.copyfile(r1_file, r1_out)
            else:
                r1_out = None
            if r2_file:
                r2_out = os.path.join(opts.workdir, os.path.basename(r2_file))
                shutil.copyfile(r2_file, r2_out)
            else:
                r2_out = None
            if i_file:
                i_out = os.path.join(opts.workdir, os.path.basename(i_file))
                shutil.copyfile(i_file, i_out)
            else:
                i_out = None
            chaste_files.append((r1_out, i_out, r2_out))

    # De-multiplex using index, keeping read1 files
    demux_files = list()
    if opts.demux_i:
        barcodes = read_barcodes(opts.barcodefile)
        print "Expected barcodes: " + ", ".join(barcodes)
        for (r1_file, i_file, r2_file) in chaste_files:
            (out_prefix_r1, out_extn) = os.path.splitext(r1_file)
            cmd = ["python", os.path.join(scriptpath, "demux.py"), "--seqs", opts.barcodefile, "--end1", r1_file, "--index", i_file]
            if r2_file:
                cmd.extend(["end2", r2_file])
            common.run_cmd(cmd)
            for seq in barcodes:
                r1_file = out_prefix_r1 + "_" + seq + out_extn
                demux_files.append((r1_file, None, None))
    else:
        demux_files = chaste_files

    # Hash and count
    if opts.verify_i:
        hash_logs = list()
        for (r1_ch_file, ind_ch_file, r2_ch_file) in chaste_files:
            parts = os.path.splitext(ind_ch_file)
            outfile = parts[0] + common.HASH_EXTENSION
            common.run_cmd(["python", os.path.join(scriptpath, "hash_index_reads.py"), "--infile", ind_ch_file, "--outfile", outfile, "--tn_end_length", str(len(opts.tn_end_seq))])
            hash_logs.append((r1_ch_file, ind_ch_file, r2_ch_file, outfile))

    # Filter first-end reads if passed hash
    if opts.verify_i:
        filtered_files = list()
        for (r1_ch_file, ind_ch_file, r2_ch_file, hash_log) in hash_logs:
            (index_ok, ratio) = check_hash(hash_log, opts.tn_end_seq)
            if index_ok:
                print "Index sequences are primarily " + opts.tn_end_seq + " (" + str(ratio * 100) + "%)"
                outfile1 = common.add_suffix(r1_ch_file, common.TNEND_SUFFIX)
                cmd = ["python", os.path.join(scriptpath, "tnend_filter.py"), "--end1", r1_ch_file, "--index", ind_ch_file, "--outfile1", outfile1, "--tn_end", opts.tn_end_seq]
                if r2_file:
                    outfile2 = common.add_suffix(r2_ch_file, common.TNEND_SUFFIX)
                    cmd.extend(["--end2", r2_ch_file, "--outfile2", outfile2])
                else:
                    outfile2 = None
                common.run_cmd(cmd)
                filtered_files.append((outfile1, None, outfile2))
            else:
                print "WARNING: most index counts do not match expected sequence - skipping Tn end filter (" + str(ratio * 100) + "% " + opts.tn_end_seq + ")" 
                filtered_files.append((r1_ch_file, None, r2_ch_file))
    else:
        filtered_files = demux_files

    # Filter and trim first-end reads
    if opts.verify_r1:
        trimmed_files = list()
        for (r1_file, i_file, r2_file) in filtered_files:
            outfile = common.add_suffix(r1_file, common.TRIM_SUFFIX)
            common.run_cmd(["python", os.path.join(scriptpath, "r1_filter.py"), "--end1", r1_file, "--outfile", outfile, "--seq", opts.tn_end_seq])
            trimmed_files.append((outfile, None, None))
    else:
        trimmed_files = filtered_files

    # De-multiplex using read2, keeping read1 and index files
    mappable_files = list()
    if opts.demux_r2:
        barcodes = read_barcodes(opts.barcodefile)
        for (r1_file, i_file, r2_file) in trimmed_files:
            (out_prefix_r1, out_extn) = os.path.splitext(r1_file)
            cmd = ["python", os.path.join(scriptpath, "demux.py"), "--seqs", opts.barcodefile, "--end1", r1_file, "--end2", r2_file]
            if i_file:
                cmd.extend(["--index", i_file])
                (out_prefix_i, _) = os.path.splitext(i_file)
            common.run_cmd(cmd)
            for seq in barcodes:
                r1_out = out_prefix_r1 + "_" + seq + out_extn
                if i_file:
                    i_out = out_prefix_i + "_" + seq + out_extn
                else:
                    i_out = None
                mappable_files.append((r1_out, i_out))
    else:
        mappable_files = [(r1, ind) for (r1, ind, r2) in trimmed_files]

    # Map primary reads against genome using desired aligner
    cores = multiprocessing.cpu_count()
    if cores > 2:
        cores -= 1
    sam_files = list()
    if opts.use_bowtie:
        print "Aligning reads with Bowtie"
        common.run_cmd([common.BOWTIE_BUILD, "--quiet", opts.reference_fa, opts.reference_fa])
        for (e1_filtered_fq, _) in mappable_files:
            parts = os.path.splitext(e1_filtered_fq)
            label = parts[0]
            sam_file = label + common.SAM_EXTENSION
            common.run_cmd([common.BOWTIE, "--sam", "--threads", str(cores),  opts.reference_fa, e1_filtered_fq, sam_file])
            sam_files.append(sam_file)
    else:
        print "Aligning reads with BWA"
        common.run_cmd([common.BWA, "index", opts.reference_fa])
        for (e1_filtered_fq, _) in mappable_files:
            parts = os.path.splitext(e1_filtered_fq)
            label = parts[0]
            aln_file = label + common.ALN_EXTENSION
            with open(aln_file, "w") as aln_fh:
#                common.run_cmd_file_out([common.BWA, "aln", "-t", str(cores), opts.reference_fa, e1_filtered_fq], aln_fh)
                common.run_cmd_file_out([common.BWA, "aln", "-l", "1000", "-t", str(cores), "-n", str(common.BWA_PCT_MISSING), opts.reference_fa, e1_filtered_fq], aln_fh)
            sam_file = label + common.SAM_EXTENSION
            with open(sam_file, "w") as sam_fh:
                common.run_cmd_file_out([common.BWA, "samse", opts.reference_fa, aln_file, e1_filtered_fq], sam_fh)
            sam_files.append(sam_file)

    # Summarize mappings
    sum_files = list()
    for sam_file in sam_files:
        parts = os.path.splitext(sam_file)
        outfile = parts[0] + common.SUM_EXTENSION
        common.run_cmd(["python", os.path.join(scriptpath, "summarize_mappings.py"), "--infile", sam_file, "--outfile", outfile])
        sum_files.append(outfile)

    # Merge slipped reads
    sum_mg_files = list()
    for sum_file in sum_files:
        outfile = common.add_suffix(sum_file, common.MERGE_SUFFIX)
        if opts.merge_slipped:
            common.run_cmd(["python", os.path.join(scriptpath, "merge_slipped.py"), "--infile", sum_file, "--outfile", outfile])
            sum_mg_files.append(outfile)
        else:
            sum_mg_files.append(sum_file)

    # Normalize read counts
    norm_files = list()
    filenum = 0
    for mg_file in sum_mg_files:
        filenum += 1
        if opts.norm_factor > 0:
            parts = os.path.splitext(mg_file)
            norm_file = parts[0] + common.NORM_SUFFIX + parts[1]
            if opts.norm_factor is None:
                common.run_cmd(["python", os.path.join(scriptpath, "norm.py"), "--infile", mg_file, "--outfile", norm_file])
            else:
                common.run_cmd(["python", os.path.join(scriptpath, "norm.py"), "--infile", mg_file, "--outfile", norm_file, "--norm-factor", str(opts.norm_factor)])
                norm_files.append(norm_file)
        else:
            norm_files.append(mg_file)

    # This file list should be pasted into the annotation command
    print "Finished reads lists: " + " ".join(norm_files)

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()

    # Input validation
    if len(args) == 0:
        print "Please supply input fastq files (must be in pairs/triplets if using index reads and/or 2nd-end reads)"
        exit(1)
    infiles = args
    for infile in infiles:
        if not os.path.isfile(infile):
            print "Input fastq " + infile + " not found"
            exit(1)
    if (opts.demux_i or opts.demux_r2) and opts.barcodefile is None:
        print "Please supply a file of expected barcodes for de-multiplexing using the -d or -e options"
        exit(1)
    if (opts.verify_r1 or opts.verify_i) and opts.tn_end_seq is None:
        print "Please provide an expected tn end sequence for verifying reads"
        exit(1)
    if opts.verify_r1 and opts.verify_i:
        print "Please choose to verify by index or verify by 1st-end read (not both)"
        exit(1)
    if opts.demux_i and opts.demux_r2:
        print "Please choose to demultiplex by index or demultiplex by 2nd-end read (not both)"
        exit(1)
    if (opts.demux_i or opts.demux_r2) and opts.barcodefile is None:
        print "Please supply a file of expected barcodes for de-multiplexing using the -d or -e options"
        exit(1)

    if not headers_ok(opts.reference_fa):
        print "WARNING: fasta headers are complex and may not be parsed correctly. We recommend simple headers: '>accession_of_genomic_element'"

    process(infiles, opts)

if __name__ == "__main__":
    main()
