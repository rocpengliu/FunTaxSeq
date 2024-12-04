#include <stdio.h>
#include <time.h>
#include <sstream>
#include "util.h"
#include "options.h"
#include "processor.h"
#include "evaluator.h"
#include "fastqreader.h"
#include "unittest.h"
#include "cmdline.h"

string command;
mutex logmtx;
int main(int argc, char* argv[]) {
    // display version info if no argument is given
    if (argc == 1) {
        cerr << "FunTaxSeq: an ultra-fast and comprehensive ." << endl << "version " << FUNTAXSEQ_VER << endl;
    }
    if (argc == 2 && strcmp(argv[1], "test") == 0) {
        UnitTest tester;
        tester.run();
        return 0;
    }
    if (argc == 2 && (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0)) {
        cerr << "funtaxseq " << FUNTAXSEQ_VER << endl;
        return 0;
    }

    cmdline::parser cmd;
    cmd.add<string>("in1", 'i', "read1 input file name", false, "");
    cmd.add<string>("in2", 'I', "read2 input file name", false, "");
    cmd.add<string>("prefix", 'X', "prefix name for output files, eg: sample01", false, "");
    cmd.add("outFReads", 0, "If specified, off-target reads will be outputed in a file");
    cmd.add<string>("out1", 'o', "file name to store read1 with on-target sequences", false, "");
    cmd.add<string>("out2", 'O', "file name to store read2 with on-target sequences", false, "");
    cmd.add<string>("samtable", 0, "sample table", false, "");
    cmd.add<string>("outdir", 0, "output directory", false, "");
    // translated search
    cmd.add<string>("tfmi", 'd', "fmi index of Protein database", false, "");
    cmd.add<string>("tmode", 'K', "searching mode either tGREEDY or tMEM (maximum exactly match). By default greedy", false, "tGREEDY");
    cmd.add<int>("mismatch", 'E', "number of mismatched amino acid in sequence comparison with protein database with default value 3", false, 3);
    cmd.add<int>("minscore", 'j', "minimum matching score of amino acid sequence in comparison with protein database with default value 65", false, 65);
    cmd.add<int>("minlength", 'J', "minimum matching length of amino acid sequence in comparison with protein database with default value 13 for GREEDY and 11 for MEM model", false, 0);
    cmd.add<int>("maxtranslength", 'm', "maximum cutoff of translated peptides, it must be no less than minlength, with default 60", false, 60);
    cmd.add<string>("codontable", 0, "select the codon table (same as blastx in NCBI), we provide 20 codon tables from 'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG31'. By default is the codontable1 (Standard Code)", false, "codontable1");
    cmd.add<string>("dbDir", 0, "dir for internal database such as ko_fullname.txt", false, "");
    // DNA search
    cmd.add<string>("dfmi", 0, "fmi index of DNA database", false, "");
    cmd.add("debug", 0, "If specified, print debug");
    // reporting
    //cmd.add<string>("json", 0, "the json format report file name", false, "fts.json");
    //cmd.add<string>("html", 'h', "the html format report file name", false, "fts.html");
    //cmd.add<string>("report_title", 'R', "should be quoted with \' or \", default is \"fts report\"", false, "fts report");
    // threading
    cmd.add<int>("thread", 'w', "worker thread number, default is 4", false, 4);
    // qother I/O
    cmd.add("phred64", '6', "indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)");
    cmd.add<int>("compression", 'z', "compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4.", false, 4);
    cmd.add("stdin", 0, "input from STDIN. If the STDIN is interleaved paired-end FASTQ, please also add --interleaved_in.");
    cmd.add("stdout", 0, "stream passing-filters reads to STDOUT. This option will result in interleaved FASTQ output for paired-end output. Disabled by default.");
    cmd.add("interleaved_in", 0, "indicate that <in1> is an interleaved FASTQ which contains both read1 and read2. Disabled by default.");
    cmd.add<int>("reads_to_process", 0, "specify how many reads/pairs to be processed. Default 0 means process all reads.", false, 0);
    cmd.add("dont_overwrite", 0, "don't overwrite existing files. Overwritting is allowed by default.");
    cmd.add("verbose", 'V', "output verbose log information (i.e. when every 1M reads are processed).");
    // adapter
    cmd.add("disable_adapter_trimming", 'A', "adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled");
    cmd.add<string>("adapter_sequence", 'a', "the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped.", false, "auto");
    cmd.add<string>("adapter_sequence_r2", 0, "the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence>", false, "auto");
    cmd.add<string>("adapter_fasta", 0, "specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file", false, "");
    cmd.add("detect_adapter_for_pe", 0, "by default, the auto-detection for adapter is for SE data input only, turn on this option to enable it for PE data.");
    // trimming
    cmd.add<int>("trim_front1", 'f', "trimming how many bases in front for read1, default is 0", false, 0);
    cmd.add<int>("trim_tail1", 't', "trimming how many bases in tail for read1, default is 0", false, 0);
    cmd.add<int>("max_len1", 'b', "if read1 is longer than max_len1, then trim read1 at its tail to make it as long as max_len1. Default 0 means no limitation", false, 0);
    cmd.add<int>("trim_front2", 'F', "trimming how many bases in front for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("trim_tail2", 'T', "trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings", false, 0);
    cmd.add<int>("max_len2", 'B', "if read2 is longer than max_len2, then trim read2 at its tail to make it as long as max_len2. Default 0 means no limitation. If it's not specified, it will follow read1's settings", false, 0);
    // polyG tail trimming
    cmd.add<int>("poly_g_min_len", 0, "the minimum length to detect polyG in the read tail. 10 by default.", false, 10);
    cmd.add("disable_trim_poly_g", 'G', "disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data");
    // polyX tail trimming
    cmd.add("trim_poly_x", 'x', "enable polyX trimming in 3' ends.");
    cmd.add<int>("poly_x_min_len", 0, "the minimum length to detect polyX in the read tail. 10 by default.", false, 10);
    // cutting by quality
    cmd.add("cut_front", '5', "move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_tail", '3', "move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.");
    cmd.add("cut_right", 'r', "move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, and then stop.");
    cmd.add<int>("cut_window_size", 'W', "the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4", false, 4);
    cmd.add<int>("cut_mean_quality", 'M', "the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20)", false, 20);
    cmd.add<int>("cut_front_window_size", 0, "the window size option of cut_front, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_front_mean_quality", 0, "the mean quality requirement option for cut_front, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_tail_window_size", 0, "the window size option of cut_tail, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_tail_mean_quality", 0, "the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified", false, 20);
    cmd.add<int>("cut_right_window_size", 0, "the window size option of cut_right, default to cut_window_size if not specified", false, 4);
    cmd.add<int>("cut_right_mean_quality", 0, "the mean quality requirement option for cut_right, default to cut_mean_quality if not specified", false, 20);
    // quality filtering
    cmd.add("disable_quality_filtering", 'Q', "quality filtering is enabled by default. If this option is specified, quality filtering is disabled");
    cmd.add<int>("qualified_quality_phred", 'q', "the quality value that a base is qualified. Default 20 means phred quality >=Q20 is qualified.", false, 20);
    cmd.add<int>("unqualified_percent_limit", 'u', "how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40%", false, 40);
    cmd.add<int>("n_base_limit", 'n', "if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5", false, 5);
    cmd.add<int>("average_qual", 'e', "if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement", false, 0);
    // length filtering
    cmd.add("disable_length_filtering", 'L', "length filtering is enabled by default. If this option is specified, length filtering is disabled");
    cmd.add<int>("length_required", 'l', "reads shorter than length_required will be discarded, default is 30.", false, 30);
    cmd.add<int>("length_limit", 0, "reads longer than length_limit will be discarded, default 0 means no limitation.", false, 0);
    // low complexity filtering
    cmd.add("low_complexity_filter", 'y', "enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).");
    cmd.add<int>("complexity_threshold", 'Y', "the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required.", false, 30);
    // filter by indexes
    cmd.add<string>("filter_by_index1", 0, "specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line", false, "");
    cmd.add<string>("filter_by_index2", 0, "specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line", false, "");
    cmd.add<int>("filter_by_index_threshold", 0, "the allowed difference of index barcode for index filtering, default 0 means completely identical.", false, 0);
    // base correction in overlapped regions of paired end data
    cmd.add("no_correction", 'C', "disable base correction in overlapped regions (only for PE data), default is enabled");
    cmd.add<int>("overlap_len_require", 0, "the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default.", false, 30);
    cmd.add<int>("overlap_diff_limit", 0, "the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default.", false, 5);
    cmd.add<int>("overlap_diff_percent_limit", 0, "the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%.", false, 20);
    // umi
    cmd.add("umi", 'U', "enable unique molecular identifier (UMI) preprocessing");
    cmd.add<string>("umi_loc", 0, "specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none", false, "");
    cmd.add<int>("umi_len", 0, "if the UMI is in read1/read2, its length should be provided", false, 0);
    cmd.add<string>("umi_prefix", 0, "if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default", false, "");
    cmd.add<int>("umi_skip", 0, "if the UMI is in read1/read2, funtaxseq can skip several bases following UMI, default is 0", false, 0);
    cmd.parse_check(argc, argv);
    if (argc == 1) {
        cerr << cmd.usage() << endl;
        return 0;
    }
    time_t t11 = time(NULL);
    stringstream ss;
    for (int i = 0; i < argc; i++) {
        ss << argv[i] << " ";
    }
    command = ss.str();
    Options * opt = new Options();
    opt->funtaxseqProgPath = string(argv[0]);
    opt->funtaxseqDir = removeStr(opt->funtaxseqProgPath, "bin/funtaxseq");
    opt->internalDBDir = cmd.get<string>("dbDir") == "" ? opt->funtaxseqDir + "database" : cmd.get<string>("dbDir");
    opt->internalDBDir = checkDirEnd(opt->internalDBDir);
    // threading
    opt->thread = cmd.get<int>("thread");
    opt->compression = cmd.get<int>("compression");
    opt->readsToProcess = cmd.get<int>("reads_to_process");
    opt->phred64 = cmd.exist("phred64");
    opt->dontOverwrite = cmd.exist("dont_overwrite");
    opt->inputFromSTDIN = cmd.exist("stdin");
    opt->outputToSTDOUT = cmd.exist("stdout");
    opt->interleavedInput = cmd.exist("interleaved_in");
    opt->verbose = cmd.exist("verbose");
    opt->debug = cmd.exist("debug");
    // adapter cutting
    opt->adapter.enabled = !cmd.exist("disable_adapter_trimming");
    opt->adapter.detectAdapterForPE = cmd.exist("detect_adapter_for_pe");
    opt->adapter.sequence = cmd.get<string>("adapter_sequence");
    opt->adapter.sequenceR2 = cmd.get<string>("adapter_sequence_r2");
    opt->adapter.fastaFile = cmd.get<string>("adapter_fasta");
    if (opt->adapter.sequenceR2 == "auto" && !opt->adapter.detectAdapterForPE && opt->adapter.sequence != "auto") {
        opt->adapter.sequenceR2 = opt->adapter.sequence;
    }
    if (!opt->adapter.fastaFile.empty()) {
        opt->loadFastaAdapters();
    }
    // trimming
    opt->trim.front1 = cmd.get<int>("trim_front1");
    opt->trim.tail1 = cmd.get<int>("trim_tail1");
    opt->trim.maxLen1 = cmd.get<int>("max_len1");
    // read2 settings follows read1 if it's not specified
    if (cmd.exist("trim_front2"))
        opt->trim.front2 = cmd.get<int>("trim_front2");
    else
        opt->trim.front2 = opt->trim.front1;

    if (cmd.exist("trim_tail2"))
        opt->trim.tail2 = cmd.get<int>("trim_tail2");
    else
        opt->trim.tail2 = opt->trim.tail1;

    if (cmd.exist("max_len2"))
        opt->trim.maxLen2 = cmd.get<int>("max_len2");
    else
        opt->trim.maxLen2 = opt->trim.maxLen1;

    // polyG tail trimming
    if (cmd.exist("disable_trim_poly_g")) {
        opt->polyGTrim.enabled = false;
    }
    opt->polyGTrim.minLen = cmd.get<int>("poly_g_min_len");

    // polyX tail trimming
    if (cmd.exist("trim_poly_x")) {
        opt->polyXTrim.enabled = true;
    }
    opt->polyXTrim.minLen = cmd.get<int>("poly_x_min_len");
    // sliding window cutting by quality
    opt->qualityCut.enabledFront = cmd.exist("cut_front");
    opt->qualityCut.enabledTail = cmd.exist("cut_tail");
    opt->qualityCut.enabledRight = cmd.exist("cut_right");
    opt->qualityCut.windowSizeShared = cmd.get<int>("cut_window_size");
    opt->qualityCut.qualityShared = cmd.get<int>("cut_mean_quality");
    if (cmd.exist("cut_front_window_size"))
        opt->qualityCut.windowSizeFront = cmd.get<int>("cut_front_window_size");
    else
        opt->qualityCut.windowSizeFront = opt->qualityCut.windowSizeShared;
    if (cmd.exist("cut_front_mean_quality"))
        opt->qualityCut.qualityFront = cmd.get<int>("cut_front_mean_quality");
    else
        opt->qualityCut.qualityFront = opt->qualityCut.qualityShared;

    if (cmd.exist("cut_tail_window_size"))
        opt->qualityCut.windowSizeTail = cmd.get<int>("cut_tail_window_size");
    else
        opt->qualityCut.windowSizeTail = opt->qualityCut.windowSizeShared;
    if (cmd.exist("cut_tail_mean_quality"))
        opt->qualityCut.qualityTail = cmd.get<int>("cut_tail_mean_quality");
    else
        opt->qualityCut.qualityTail = opt->qualityCut.qualityShared;

    if (cmd.exist("cut_right_window_size"))
        opt->qualityCut.windowSizeRight = cmd.get<int>("cut_right_window_size");
    else
        opt->qualityCut.windowSizeRight = opt->qualityCut.windowSizeShared;
    if (cmd.exist("cut_right_mean_quality"))
        opt->qualityCut.qualityRight = cmd.get<int>("cut_right_mean_quality");
    else
        opt->qualityCut.qualityRight = opt->qualityCut.qualityShared;

    // raise a warning if cutting option is not enabled but -W/-M is enabled
    if (!opt->qualityCut.enabledFront && !opt->qualityCut.enabledTail && !opt->qualityCut.enabledRight) {
        if (cmd.exist("cut_window_size") || cmd.exist("cut_mean_quality")
                || cmd.exist("cut_front_window_size") || cmd.exist("cut_front_mean_quality")
                || cmd.exist("cut_tail_window_size") || cmd.exist("cut_tail_mean_quality")
                || cmd.exist("cut_right_window_size") || cmd.exist("cut_right_mean_quality"))
            cerr << "WARNING: you specified the options for cutting by quality, but forogt to enable any of cut_front/cut_tail/cut_right. This will have no effect." << endl;
    }
    // quality filtering
    opt->qualfilter.enabled = !cmd.exist("disable_quality_filtering");
    opt->qualfilter.qualifiedQual = num2qual(cmd.get<int>("qualified_quality_phred"));
    opt->qualfilter.unqualifiedPercentLimit = cmd.get<int>("unqualified_percent_limit");
    opt->qualfilter.avgQualReq = cmd.get<int>("average_qual");
    opt->qualfilter.nBaseLimit = cmd.get<int>("n_base_limit");
    // length filtering
    opt->lengthFilter.enabled = !cmd.exist("disable_length_filtering");
    opt->lengthFilter.requiredLength = cmd.get<int>("length_required");
    opt->lengthFilter.maxLength = cmd.get<int>("length_limit");
    //get trans mode first;
    opt->mTransSearchOptions->tmode = cmd.get<string>("tmode");
    if (cmd.get<int>("minlength") == 0) {
        if (opt->mTransSearchOptions->mode == tGREEDY) {
            opt->mTransSearchOptions->minAAFragLength = 13;
        } else {
            opt->mTransSearchOptions->minAAFragLength = 10;
        }
    } else {
        opt->mTransSearchOptions->minAAFragLength = cmd.get<int>("minlength");
    }
    opt->lengthFilter.requiredLength = max(opt->lengthFilter.requiredLength, static_cast<int> (opt->mTransSearchOptions->minAAFragLength) * 3);
    opt->lengthFilter.maxLength = cmd.get<int>("length_limit");
    // low complexity filter
    opt->complexityFilter.enabled = cmd.exist("low_complexity_filter");
    opt->complexityFilter.threshold = (min(100, max(0, cmd.get<int>("complexity_threshold")))) / 100.0;
    // overlap correction
    opt->correction.enabled = !cmd.exist("no_correction");
    opt->overlapRequire = cmd.get<int>("overlap_len_require");
    opt->overlapDiffLimit = cmd.get<int>("overlap_diff_limit");
    opt->overlapDiffPercentLimit = cmd.get<int>("overlap_diff_percent_limit");
    // umi
    opt->umi.enabled = cmd.exist("umi");
    opt->umi.length = cmd.get<int>("umi_len");
    opt->umi.prefix = cmd.get<string>("umi_prefix");
    opt->umi.skip = cmd.get<int>("umi_skip");
    if (opt->umi.enabled) {
        string umiLoc = cmd.get<string>("umi_loc");
        str2lower(umiLoc);
        if (umiLoc.empty())
            error_exit("You've enabled UMI by (--umi), you should specify the UMI location by (--umi_loc)");
        if (umiLoc != "index1" && umiLoc != "index2" && umiLoc != "read1" && umiLoc != "read2" && umiLoc != "per_index" && umiLoc != "per_read") {
            error_exit("UMI location can only be index1/index2/read1/read2/per_index/per_read");
        }

        //if (!opt->isPaired() && (umiLoc == "index2" || umiLoc == "read2"))
            //error_exit("You specified the UMI location as " + umiLoc + ", but the input data is not paired end.");
        if (opt->umi.length == 0 && (umiLoc == "read1" || umiLoc == "read2" || umiLoc == "per_read"))
            error_exit("You specified the UMI location as " + umiLoc + ", but the length is not specified (--umi_len).");
        if (umiLoc == "index1") {
            opt->umi.location = UMI_LOC_INDEX1;
        } else if (umiLoc == "index2") {
            opt->umi.location = UMI_LOC_INDEX2;
        } else if (umiLoc == "read1") {
            opt->umi.location = UMI_LOC_READ1;
        } else if (umiLoc == "read2") {
            opt->umi.location = UMI_LOC_READ2;
        } else if (umiLoc == "per_index") {
            opt->umi.location = UMI_LOC_PER_INDEX;
        } else if (umiLoc == "per_read") {
            opt->umi.location = UMI_LOC_PER_READ;
        }
    }
    opt->mEvaluation.supportEvaluation = (!opt->inputFromSTDIN && opt->in1 != "/dev/stdin");
    opt->mEvaluation.disable_trim_poly_g = cmd.exist("disable_trim_poly_g");

    opt->outdir = cmd.get<string>("outdir");
    opt->sampleTable = cmd.get<string>("samtable");
    if(opt->sampleTable.empty()) {
        Sample sam;
        sam.ff = cmd.get<string>("in1");
        sam.rf = cmd.get<string>("in2");
        sam.prefix = cmd.get<string>("prefix");
        opt->outFRFile = sam.prefix + "_funtax.txt.gz";
        opt->jsonFile = sam.prefix + ".json";
        opt->htmlFile = sam.prefix + ".html";
        opt->reportTitle = sam.prefix;
        opt->samVec.push_back(sam);
    } else {
        opt->parseSampleTable();
    }

    //translated search
    opt->mTransSearchOptions->tCodonTable = cmd.get<string>("codontable");
    opt->deterCodonTable();

    opt->mTransSearchOptions->misMatches = cmd.get<int>("mismatch");
    opt->mTransSearchOptions->minScore = cmd.get<int>("minscore");
    opt->mTransSearchOptions->maxTransLength = cmd.get<int>("maxtranslength");
    opt->mTransSearchOptions->maxTransLength = max(opt->mTransSearchOptions->maxTransLength, opt->mTransSearchOptions->minAAFragLength);
    opt->mTransSearchOptions->maxTransLength = min((unsigned) 60, opt->mTransSearchOptions->maxTransLength);
    opt->mTransSearchOptions->tfmi = cmd.get<string>("tfmi");
    opt->mDNASearchOptions->dfmi = cmd.get<string>("dfmi");
    BwtFmiDBPair* mBwtfmiDBPair = new BwtFmiDBPair(opt);
     for(int i = 0; i < opt->samVec.size(); ++i){
         time_t t1 = time(NULL);
         opt->mHomoSearchOptions->reset2Default();
         auto sam = opt->samVec.at(i);
         cerr << "Start to process sample " << sam.prefix << ", " << (i + 1) << " out of " << opt->samVec.size() << " samples" << endl;
         opt->in1 = sam.ff;
         opt->in2 = sam.rf;
         opt->prefix = sam.prefix;
         opt->outFRFile = opt->prefix + "_funtax.txt.gz";
         opt->jsonFile = opt->prefix + ".json";
         opt->htmlFile = opt->prefix + ".html";
         opt->reportTitle = opt->prefix;
         if (opt->debug)
             cCout(opt->prefix + " " + opt->in1 + " " + opt->in2);
         Evaluator eva(opt);
         if (opt->mEvaluation.supportEvaluation){
             eva.evaluateSeqLen();
         }
         long readNum = 0;
         if (opt->shallDetectAdapter(false)){
             if (!opt->mEvaluation.supportEvaluation){
                 cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
             } else{
                 cerr << "Detecting adapter sequence for read1..." << endl;
                 string adaptA = eva.evalAdapterAndReadNum(readNum, false);
                 if (adaptA.length() > 60)
                     adaptA.resize(0, 60);
                 if (adaptA.length() > 0){
                     opt->adapter.sequence = adaptA;
                     opt->adapter.detectedAdapter1 = adaptA;
                 } else {
                     cerr << "No adapter detected for read1" << endl;
                     opt->adapter.sequence = "";
                 }
                 cerr << endl;
             }
         }
         if (opt->shallDetectAdapter(true)){
             if (!opt->mEvaluation.supportEvaluation){
                 cerr << "Adapter auto-detection is disabled for STDIN mode" << endl;
             } else{
                 cerr << "Detecting adapter sequence for read2..." << endl;
                 string adaptA = eva.evalAdapterAndReadNum(readNum, true);
                 if (adaptA.length() > 60)
                     adaptA.resize(0, 60);
                 if (adaptA.length() > 0){
                     opt->adapter.sequenceR2 = adaptA;
                     opt->adapter.detectedAdapter2 = adaptA;
                 } else {
                     cerr << "No adapter detected for read2" << endl;
                     opt->adapter.sequenceR2 = "";
                 }
                 cerr << endl;
             }
         }

        if (!opt->mEvaluation.disable_trim_poly_g && opt->mEvaluation.supportEvaluation){
             bool twoColorSystem = eva.isTwoColorSystem();
             if (twoColorSystem){
                 opt->polyGTrim.enabled = true;
             }
        }
        Processor* p = new Processor(opt);
        p->process(mBwtfmiDBPair);
        if (p) {
            delete p;
            p = NULL;
        }
        opt->mHomoSearchOptions->calculate();
        time_t t2 = time(NULL);
        cerr << endl << "JSON report: " << opt->jsonFile << endl;
        cerr << "HTML report: " << opt->htmlFile << endl;
        cerr << "funtaxseq v" << FUNTAXSEQ_VER << ", time used: " << convertSeconds((t2) - t1) 
        << " mapped " << opt->mHomoSearchOptions->mappedReads << " reads(" 
        << opt->mHomoSearchOptions->mappedReadPer << "%) for sample: " << sam.prefix << "!" << endl;
    }
    if(mBwtfmiDBPair){
        delete mBwtfmiDBPair;
        mBwtfmiDBPair = nullptr;
    }
    time_t t22 = time(NULL);
    cerr << endl << command << endl;
    cerr << "funtaxdecoder v" << FUNTAXSEQ_VER << " time used: " << convertSeconds((t22) - t11) << " processed " << opt->samVec.size() << " samples" << endl;
    if(opt) delete opt; opt = nullptr;
    return 0;
}