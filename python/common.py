#!/usr/bin/env python
# Common constants and methods for the Tn-seq pipeline
#
# Copyright (c) 2014 University of Washington

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
import os.path
import subprocess

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
IDX_RATIO = .5
BWA = "/usr/bin/bwa"
BWA_SEED_DIFF = 2
BWA_PCT_MISSING = .035
TN_END_LENGTH = 6
BOWTIE = "/usr/bin/bowtie"
BOWTIE_BUILD = "/usr/bin/bowtie-build"
ALN_EXTENSION = ".sai"
SAM_EXTENSION = ".sam"
SUM_EXTENSION = "_sum.txt"
HASH_EXTENSION = ".index.log"
CHASTE_SUFFIX = "_ch"
TNEND_SUFFIX = "_iPass"
TRIM_SUFFIX = "_trim"
MERGE_SUFFIX = "_mg"
NORM_SUFFIX = "_norm"
NORM_FACTOR = 10000000
WORKING_DIR = "work"

ALL_SUFFIX = "_all"
Q0_SUFFIX = "_q0"
ANNO_SUFFIX = "_Annot"
TAB_SUFFIX = "_HitsToTab"
READSCOMP = "reads_cmp.txt"
OUTFILE_ANNO = "AnnotatedHits.txt"
OUTFILE_TAB = "HitsPerGene.txt"

#------------------------------------------------------------------------------
# methods
#------------------------------------------------------------------------------
# Add a suffix to a filename before the file extension
def add_suffix(filename, suffix):
    parts = os.path.splitext(filename)
    new_name = parts[0] + suffix + parts[1]
    return new_name

# Check and run a command or exit with an error
def run_cmd(cmd):
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        badcmd = " ".join(e.cmd)
        print "Error running command: " + badcmd
        exit(1)

# Check and run a command or exit with an error, piping stdout to the specified file
def run_cmd_file_out(cmd, stdout_fh):
    try:
        subprocess.check_call(cmd, stdout=stdout_fh)
    except subprocess.CalledProcessError as e:
        badcmd = " ".join(e.cmd)
        print "Error running command: " + badcmd
        exit(1)
