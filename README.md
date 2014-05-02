Manual for Tn-seq - An analsys pipeline for the Tn-seq protocol

Overview
---------------------
The Tn-seq pipeline is run using two master scripts. The first step generates
lists of reads per location and can be run seperately for multiple Tn-seq runs.
The second step annotates the locations and tabulates hits per gene and can
incorporate input files from multiple Tn-seq runs for comparison.

Software requirements
---------------------
Python - we use version 2.6.6 but any 2.x version above that should work.

Read mapping software:
BWA (tested with version 0.7.4) or Bowtie (tested with version 0.12.7)

Running the scripts
---------------------
```
Step 1: mapping (process_map.py)
Usage: process_map.py [options] firstend_fastq_1 index_fastq_1 secondend_fastq_1 ...
Options:
  -h, --help                show this help message and exit
  -r, --reffile             path to reference genome fasta
  -j, --tn_verify_by_read1  use 1st-end read to verify transposon end (default: False)
  -i, --tn_verify_by_index  use index read to verify transposon end (default: False)
  -t, --tn_end              expected transposon end sequence
  -d, --demux_index         use index read to demultiplex (default: False)
  -e, --demux_read2         use 2nd end read to demultiplex (default: False)
  -b, --barcodefile         path to file listing expected barcode sequences
  -c, --chastity            run chastity filter (default: False)
  -n, --normfactor          read count normalization factor (default: 10,000,000)
                            (0 = don't normalize)
  -s, --merge_slipped       merge slipped reads (default: False)
  -u, --use_bowtie          map reads using Bowtie (default: use BWA)
  -w, --workingdir          working directory for input and output files (default: work)

Step 2: annotating (process_annotate_tablulate.py)
Usage: process_annotate_tabulate.py [options] reads_list_1 reads_list_2 ...
Options:
  -h, --help            show this help message and exit
  -r, --reffile         path to reference genome fasta
  -a, --annofiles       path to reference .ptt annotation file(s) (comma-separated list if
			using more than one; order must match sequences in reference fasta)
  -o, --outfile_anno    path to final annotated output file (default: work/AnnotatedHits.txt)
  -p, --outfile_tab     path to final counts tabulated by gene (default: work/HitsPerGene.txt)
  -w, --workingdir      working directory for input and output files (default: work)
```

Examples
---------------------
```
python process_map.py --barcodefile barcodes.txt --chastity --demux_read2 --tn_verify_by_index --reffile combined.fna --tn_end AGACAG --workingdir work r1.fq ind.fq r2.fq

python process_annotate_tabulate.py --annofiles CP000086.ptt,CP000085.ptt --reffile combined.fna --workingdir work work/r1_ch_iPass_ACGTGA_sum_norm.txt work/r1_ch_iPass_CTAGTG_sum_norm.txt work/r1_ch_iPass_GATCAC_sum_norm.txt work/r1_ch_iPass_TGCACT_sum_norm.txt
```

Additional notes
---------------------
The comma-separated list of .ptt annotation files should not have spaces between the
files (only a comma).
If your reference genome contains multiple replicons, combine their fasta files into
a single fasta before running this software. The order of the sequences should be
the same as the order of the comma-separated list of .ptt files.
