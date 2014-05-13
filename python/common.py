#!/usr/bin/env python
# Common constants and methods for the Tn-seq pipeline
#
# Copyright (c) 2014 University of Washington

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
import os.path
import subprocess
import re

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
IDX_RATIO = .5
BWA = "/usr/bin/bwa"
BWA_SEED_DIFF = 2
BWA_PCT_MISSING = .035
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
SCRIPTS_PATH = "."

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

# Replicon names are in the fasta headers - read fasta and return a list of names
def read_replicon_names(fasta):
    replicons = dict()
    repl_num = 0
    with open(fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                repl_long = line[1:-1]
                repl = repl_long.split(" ", 1)[0]
                replicons[repl_num] = repl
                repl_num += 1
    return replicons

# Read annotation data from one or more .ptt files and return a hash
def read_annotations(annofile_list, replicon_list):
    annotations = dict()
    count = 0
    for filenum, annofile in enumerate(annofile_list):
        replicon = replicon_list[filenum]
        annotations[replicon] = dict()
        with open(annofile, "r") as fh:
            for line in fh:
                mobj = re.match("^(\d+)\.\.(\d+)", line)
                if mobj:
                    startpos = mobj.group(1)
                    endpos = mobj.group(2)
                    (loc, strand, length, pid, gene, synonym, code, cog, product) = line.rstrip().split('\t')
                    annotations[replicon][pid] = dict();
                    annotations[replicon][pid]['locus_tag'] = synonym
                    annotations[replicon][pid]['startpos'] = int(startpos)
                    annotations[replicon][pid]['endpos'] = int(endpos)
                    annotations[replicon][pid]['strand'] = strand
                    annotations[replicon][pid]['length'] = length
                    annotations[replicon][pid]['info'] = '\t'.join([gene, code, cog, product])
                    count += 1
    print "read " + str(count) + " annotations for " + str(len(annotations.keys())) + " replicon(s)"
    return annotations
