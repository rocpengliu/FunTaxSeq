# FunTaxSeq
## a unified platform for simultaneous taxonomic and functional profilings based on tiered search from metagenomic sequence data!

### Download and build funtaxseq
```shell
# get source (you can also use browser to download from main or releases)
git clone https://github.com/rocpengliu/FunTaxSeq.git

# build
cd FunTaxSeq
make clean && make;

# after the building, there are 4 exe files in the dir bin
# fts is used to map reads
# ftd is to generate taxa and gene abundance tables
# mkbwt and mkfmi are used to build databases
```

### sample.txt separated by tab
```
out_fts/sim_rep1    raw/sim_rep1.R1_sub.fastq.gz    raw/sim_rep1.R2_sub.fastq.gz
out_fts/sim_rep2    raw/sim_rep2.R1_sub.fastq.gz    raw/sim_rep2.R2_sub.fastq.gz
out_fts/sim_rep3    raw/sim_rep3.R1_sub.fastq.gz    raw/sim_rep3.R2_sub.fastq.gz
```

## first try
```
FunTaxSeq/bin/fts --samtable sample.txt --dfmi db/dna_sub.fmi --tfmi db/pro_sub.fmi -w 12 -V
FunTaxSeq/bin/ftd --samdir out_fts -b db/ -o out_ftd/out_ftd -w 12 -V
```

### all options for fts
```
usage: ./fts [options] ...
options:
  -i, --in1                              read1 input file name (string [=])
  -I, --in2                              read2 input file name (string [=])
  -X, --prefix                           prefix name for output files, eg: sample01 (string [=])
  -o, --out1                             file name to store read1 with on-target sequences (string [=])
  -O, --out2                             file name to store read2 with on-target sequences (string [=])
      --samtable                         sample table (string [=])
      --outdir                           output directory (string [=])
  -d, --tfmi                             fmi index of Protein database (string [=])
  -K, --tmode                            searching mode either GREEDY or MEM (maximum exactly match). By default greedy (string [=GREEDY])
  -E, --tmismatch                        number of mismatched amino acid in sequence comparison with protein database with default value 3 (int [=3])
  -j, --tminscore                        minimum matching score of amino acid sequence in comparison with protein database with default value 65 (int [=65])
  -J, --tminlength                       minimum matching length of amino acid sequence in comparison with protein database with default value 13 for GREEDY and 11 for MEM model (int [=0])
      --tlenper                          percentage of matching amino acid in sequence comparison with default value 0.8 (float [=0.8])
      --codontable                       select the codon table (same as blastx in NCBI), we provide 20 codon tables from 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1'. By default is the codontable1 (Standard Code) (string [=codontable1])
      --dfmi                             fmi index of DNA database (string [=])
      --dmode                            searching mode either GREEDY or MEM (maximum exactly match) in a DNA database. By default greedy (string [=GREEDY])
      --dmismatch                        number of mismatches in sequence comparison with DNA database with default value 6 (int [=6])
      --dminlength                       minimum matching length of DNA sequence in comparison with DNA database with default value 60 (int [=60])
      --dlenper                          percentage of matching DNA in sequence comparison with default value 0.9 (float [=0.9])
      --dminscore                        minimum matching score of DNA sequence in comparison with DNA database with default value 50 (int [=50])
      --hfmi                             fmi index of host DNA database (string [=])
      --hmode                            searching mode either GREEDY or MEM (maximum exactly match) in a host database. By default greedy (string [=GREEDY])
      --hmismatch                        number of mismatches in sequence comparison with host DNA database with default value 5 (int [=5])
      --hminlength                       minimum matching length against host DNA database with default value 95 (int [=95])
      --hlenper                          percentage of matching DNA in sequence comparison with default value 0.90 (float [=0.9])
      --hminscore                        minimum matching score against host DNA database with default value 95 (int [=95])
      --mfmi                             fmi index of marker genes DNA database (string [=])
      --mmode                            searching mode either GREEDY or MEM (maximum exactly match) in marker genes (16s, ITS) DNA database. By default greedy (string [=GREEDY])
      --mmismatch                        number of mismatches in sequence comparison with marker genes DNA database with default value 3 (int [=3])
      --mminlength                       minimum matching length against marker genes DNA database with default value 95 (int [=95])
      --mlenper                          percentage of matching DNA in sequence comparison with default value 0.95 (float [=0.97])
      --mminscore                        minimum matching score against marker genes DNA database with default value 95 (int [=95])
      --debug                            If specified, print debug
  -w, --thread                           number of worker thread, default is 4 (int [=4])
  -6, --phred64                          indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
  -z, --compression                      compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4])
      --stdin                            input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.
      --stdout                           stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. Disabled by default.
      --interleaved_in                   indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.
      --reads_to_process                 specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])
      --dont_overwrite                   don't overwrite existing files. Overwritting is allowed by default.
  -V, --verbose                          output verbose log information (i.e. when every 1M reads are processed).
  -A, --disable_adapter_trimming         adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
  -a, --adapter_sequence                 the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
      --adapter_sequence_r2              the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=auto])
      --adapter_fasta                    specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
      --detect_adapter_for_pe            by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.
  -f, --trim_front1                      trimming how many bases in front for read1, default is 0 (int [=0])
  -t, --trim_tail1                       trimming how many bases in tail for read1, default is 0 (int [=0])
  -b, --max_len1                         if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation (int [=0])
  -F, --trim_front2                      trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
  -T, --trim_tail2                       trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])
  -B, --max_len2                         if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings (int [=0])
      --poly_g_min_len                   the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
  -G, --disable_trim_poly_g              disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data
  -x, --trim_poly_x                      enable polyX trimming in 3' ends.
      --poly_x_min_len                   the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
  -5, --cut_front                        move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
  -3, --cut_tail                         move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
  -r, --cut_right                        move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.
  -W, --cut_window_size                  the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
  -M, --cut_mean_quality                 the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])
      --cut_front_window_size            the window size option of cut_front, default to cut_window_size if not specified (int [=4])
      --cut_front_mean_quality           the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
      --cut_tail_window_size             the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
      --cut_tail_mean_quality            the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
      --cut_right_window_size            the window size option of cut_right, default to cut_window_size if not specified (int [=4])
      --cut_right_mean_quality           the mean quality requirement option for cut_right, default to cut_mean_quality if not specified (int [=20])
  -Q, --disable_quality_filtering        quality filtering is enabled by default. If this option is specified, quality filtering is disabled
  -q, --qualified_quality_phred          the quality value that a base is qualified. Default 20 means phred quality >=Q20 is qualified. (int [=20])
  -u, --unqualified_percent_limit        how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
  -n, --n_base_limit                     if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
  -e, --average_qual                     if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])
  -L, --disable_length_filtering         length filtering is enabled by default. If this option is specified, length filtering is disabled
  -l, --length_required                  reads shorter than length_required will be discarded, default is 30. (int [=30])
      --length_limit                     reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])
  -y, --disable_low_complexity_filter    disable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]). Disabled if specified
  -Y, --complexity_threshold             the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])
      --filter_by_index1                 specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
      --filter_by_index2                 specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
      --filter_by_index_threshold        the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])
  -C, --no_correction                    disable base correction in overlapped regions (only for PE data), default is enabled
      --overlap_len_require              the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
      --overlap_diff_limit               the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
      --overlap_diff_percent_limit       the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])
  -U, --umi                              enable unique molecular identifier (UMI) preprocessing
      --umi_loc                          specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
      --umi_len                          if the UMI is in read1/read2, its length should be provided (int [=0])
      --umi_prefix                       if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
      --umi_skip                         if the UMI is in read1/read2, funtaxseq can skip several bases following UMI, default is 0 (int [=0])
  -?, --help                             print this message
```

### all options for ftd to decode the output of FunTaxSeq to taxonomic and functional profilings.
```
usage: ./ftd [options] ...
options:
  -s, --samdir       sample directory (string [=])
  -o, --outprefix    output file prefix (string [=])
  -b, --database     database directory (string [=])
  -g, --useogtree    If specified, using ortholog tree(slow)
  -w, --thread       number of worker thread, default is 4 (int [=4])
  -d, --debug        If specified, print debug
  -V, --verbose      output verbose
  -?, --help         print this message
```

### TO DO LIST
```

```
