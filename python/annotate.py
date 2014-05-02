#!/usr/bin/env python
# Purpose: Annotate transposon insertion locations using
# one or more .ptt files for the genome
#
# Copyright (c) 2014 University of Washington

#-----------------------------------------------------
# modules
#-----------------------------------------------------
import sys
import optparse
import os
import re
import glob

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
PROGRESS_FREQ = 100
OUT_SUFFIX = "_annot.xls"

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(prog=sys.argv[0],
                          usage=usage,
                          add_help_option=True)
    parser.add_option("-i", "--infile", action="store", type="string", dest="infile", help="path to read counts summary file")
    parser.add_option("-a", "--annofiles", action="store", type="string", dest="annofiles", help="path to .ptt annotation file(s) (comma-separated list if using more than one; order must match sequences in reference fasta)")
    parser.add_option("-f", "--fasta", action="store", type="string", dest="fasta", help="path to reference fasta containing replicon names")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="path to output file with merged read counts (default <infile>" + OUT_SUFFIX + ")")

    opts, args = parser.parse_args()
    if (opts.infile is None or opts.annofiles is None or opts.fasta is None):
        parser.print_help()
        exit(1)

    # Create the default output filename if not provided
    if (opts.outfile is None):
        parts = os.path.splitext(os.path.basename(opts.infile))
        opts.outfile = parts[0] + OUT_SUFFIX + parts[1]

    return opts, args

# Replicon names are in the fasta headers
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

def read_annotations(annofile_list, replicon_list):
    annotations = dict()
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
                    annotations[replicon][pid]['info'] = '\t'.join([gene, code, cog, product])
    print "read annotation files"
    return annotations

def annotate(replicon, position, direction, annotations):
    # Positions immediately right and left of the Tn
    if direction == 'F':
        lpos = position
        rpos = position - 8
    else:
        lpos = position + 8
        rpos = position

    effectivepos = "";
    locus_tag = "";
    startpos = "";
    endpos = "";
    strand = "";
    info = "\t\t\t";
    relativepos = "";
    notes = "";

    # Find all orfs possibly interrupted
    interrupted_genes = list();
    for pid, pinfo in annotations[replicon].iteritems():
        if ((pinfo['strand'] == '+' and pinfo['startpos'] <= lpos and pinfo['endpos'] >= lpos) or
        (pinfo['strand'] == '-' and pinfo['startpos'] <= rpos and pinfo['endpos'] >= rpos)):
            interrupted_genes.append(pid)

    if len(interrupted_genes) == 1:
        # One gene affected
        pid = interrupted_genes[0]
        pinfo = annotations[replicon][pid]
        if pinfo['strand'] == '+':
            effectivepos = lpos
            relativepos = str(1 + lpos - pinfo['startpos']) + "(" + str(1 + pinfo['endpos'] - pinfo['startpos']) + ")"
        else:
            effectivepos = rpos
            relativepos = str(1 + pinfo['endpos'] - rpos) + "(" + str(1 + pinfo['endpos'] - pinfo['startpos']) + ")"
        (gene_id, locus_tag, startpos, endpos, strand, info) = (pid, pinfo['locus_tag'], pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['info'])

    elif len(interrupted_genes) == 0:
        # Intergenic
        notes = "intergenic"
        annos_by_pos = sorted(annotations[replicon].iteritems(), key=lambda (k,v): v['startpos'])
        (orf_to_left, _) = annos_by_pos[0]
        (orf_to_right, _) = annos_by_pos[-1]
        prev_pid = orf_to_left
        for anno in annos_by_pos:
            (pid, _) = anno
            if annotations[replicon][pid]['startpos'] > rpos:
                orf_to_left = prev_pid
                orf_to_right = pid
                break
            prev_pid = pid

        if annotations[replicon][orf_to_left]['strand'] == annotations[replicon][orf_to_right]['strand']:
            # Flanking orfs on same strand
            if annotations[replicon][orf_to_left]['strand'] == "+":
                pinfo = annotations[replicon][orf_to_right]
                (effectivepos, gene_id, locus_tag, startpos, endpos, strand, info) = (lpos, orf_to_right, pinfo['locus_tag'], pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['info'])
                relativepos = str(0 - (pinfo['startpos'] - lpos))
            else:
                pinfo = annotations[replicon][orf_to_left]
                (effectivepos, gene_id, locus_tag, startpos, endpos, strand, info) = (rpos, orf_to_left, pinfo['locus_tag'], pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['info'])
                relativepos = str(0 - (rpos - pinfo['endpos']))
        else:
            # Flanking orfs on opposite strands
            if rpos - annotations[replicon][orf_to_left]['endpos'] <= annotations[replicon][orf_to_right]['startpos'] - lpos:
                # Closer to gene on the left
                pinfo = annotations[replicon][orf_to_left]
                if pinfo['strand'] == "+":
                    (effectivepos, gene_id, locus_tag, startpos, endpos, strand, info) = (lpos, orf_to_left, pinfo['locus_tag'], pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['info'])
                    relativepos = str(lpos - pinfo['endpos'])
                else:
                    (effectivepos, gene_id, locus_tag, startpos, endpos, strand, info) = (rpos, orf_to_left, pinfo['locus_tag'], pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['info'])
                    relativepos = str(0 - (rpos - pinfo['endpos']))
            else:
                # Closer to gene on the right
                pinfo = annotations[replicon][orf_to_right]
                if pinfo['strand'] == "+":
                    (effectivepos, gene_id, locus_tag, startpos, endpos, strand, info) = (lpos, orf_to_right, pinfo['locus_tag'], pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['info'])
                    relativepos = (0 - (pinfo['startpos'] - lpos))
                else:
                    (effectivepos, gene_id, locus_tag, startpos, endpos, strand, info) = (rpos, orf_to_right, pinfo['locus_tag'], pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['info'])
                    relativepos = str(pinfo['startpos'] - rpos)
    elif len(interrupted_genes) == 2:
        # 2 overlapping genes - choose the one with the earliest relative interruption
        notes += "pos. 2 genes aff."
        interrupted_gene_locs = [(annotations[replicon][pid]['startpos'], pid) for pid in interrupted_genes]
        (_, orf_to_left) = min(interrupted_gene_locs, key = lambda k: k[0])
        (_, orf_to_right) = max(interrupted_gene_locs, key = lambda k: k[0])
        if annotations[replicon][orf_to_left]['strand'] == "+":
            relpos_l = lpos - annotations[replicon][orf_to_left]['startpos'] + 1
        else:
            relpos_l = annotations[replicon][orf_to_left]['endpos'] - rpos + 1
        if annotations[replicon][orf_to_right]['strand'] == "+":
            relpos_r = lpos - annotations[replicon][orf_to_right]['startpos'] + 1
        else:
            relpos_r = annotations[replicon][orf_to_right]['endpos'] - rpos + 1
        if relpos_l <= relpos_r:
            # Choose the left one
            pinfo = annotations[replicon][orf_to_left]
            if pinfo['strand'] == "+":
                effectivepos = lpos
            else:
                effectivepos = rpos
            relativepos = str(relpos_l) + "(" + str(pinfo['endpos'] - pinfo['startpos'] + 1) + ")"
            (gene_id, locus_tag, startpos, endpos, strand, info) = (orf_to_left, pinfo['locus_tag'], pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['info'])
        else:
            # Choose the right one
            pinfo = annotations[replicon][orf_to_right]
            if pinfo['strand'] == "+":
                effectivepos = lpos
            else:
                effectivepos = rpos
            relativepos = str(relpos_r) + "(" + str(pinfo['endpos'] - pinfo['startpos'] + 1) + ")"
            (gene_id, locus_tag, startpos, endpos, strand, info) = (orf_to_right, pinfo['locus_tag'], pinfo['startpos'], pinfo['endpos'], pinfo['strand'], pinfo['info'])
    else:
        gene_id = "multiple possible"
        notes = "multiple genes potentially affected"

    if effectivepos != "" and effectivepos != position:
        if notes != "":
            notes += "; "
        notes += "ef.pos.pred."

    return "\t".join(map(str, [replicon, position, direction, effectivepos, gene_id, locus_tag, startpos, endpos, strand, info, relativepos, notes])) + "\n";

def write_annotated(infile, outfile, annotations):
    count = 0
    with open(infile, "r") as infh:
        with open(outfile, "w") as outfh:
            in_header = infh.readline()
            out_header = "Replicon\tPosition of insertion(t+1)\tDirection\tEffective Position\tPID\tLocus\tFrom\tTo\tStrand\tGene\tCode\tCOG\tProduct\tPosition relative to locus\tNotes\n"
            outfh.write(out_header)
            for line in infh:
                (replicon, position, direction, other) = line.rstrip().split("\t", 3)
                outline = annotate(replicon, int(position), direction, annotations)
                outfh.write(outline)
                count += 1
                if count % PROGRESS_FREQ == 0:
                    sys.stdout.write("\rannotated " + str(count) + " positions")
                    sys.stdout.flush()
    print "\rannotated " + str(count) + " positions"
    return count

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()
    print "Annotating " + opts.infile

    all_replicons = read_replicon_names(opts.fasta)

    # Separate multiple annotation files
    annofile_list = opts.annofiles.split(",")
    if len(annofile_list) != len(all_replicons):
        print "ERROR: wrong number of replicon names for replicon files"
        exit(1)
    print "Using " + str(len(annofile_list)) + " .ptt files"

    annotations = read_annotations(annofile_list, all_replicons)

    lines_written = write_annotated(opts.infile, opts.outfile, annotations)
    print str(lines_written) + " lines written to " + opts.outfile

if __name__ == "__main__":
    main()
