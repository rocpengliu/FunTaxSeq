#ifndef OPTIONS_H
#define OPTIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string.h>
#include <functional>

#include "fastareader.h"
#include "util.h"
#include "common.h"


using namespace std;

#define UMI_LOC_NONE 0
#define UMI_LOC_INDEX1 1
#define UMI_LOC_INDEX2 2
#define UMI_LOC_READ1 3
#define UMI_LOC_READ2 4
#define UMI_LOC_PER_INDEX 5
#define UMI_LOC_PER_READ 6

class DuplicationOptions {
public:
    DuplicationOptions() {
        enabled = true;
        keylen = 12;
        histSize = 32;
    }
public:
    bool enabled;
    int keylen;
    int histSize;
};

class LowComplexityFilterOptions {
public:
    LowComplexityFilterOptions() {
        enabled = true;
        threshold = 0.3;
    }
public:
    bool enabled;
    double threshold;
};

class PolyGTrimmerOptions {
public:
    PolyGTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class PolyXTrimmerOptions {
public:
    PolyXTrimmerOptions() {
        enabled = false;
        minLen = 10;
    }
public:
    bool enabled;
    int minLen;
};

class UMIOptions {
public:
    UMIOptions() {
        enabled = false;
        location = UMI_LOC_NONE;
        length = 0;
        skip = 0;
    }
public:
    bool enabled;
    int location;
    int length;
    int skip;
    string prefix;
    string separator;
};

class CorrectionOptions {
public:
    CorrectionOptions() {
        enabled = true;
    }
public:
    bool enabled;
};

class QualityCutOptions {
public:
    QualityCutOptions() {
        enabledFront = false;
        enabledTail = false;
        enabledRight = false;
        windowSizeShared = 4;
        qualityShared = 20;
        windowSizeFront = windowSizeShared;
        qualityFront = qualityShared;
        windowSizeTail = windowSizeShared;
        qualityTail = qualityShared;
        windowSizeRight = windowSizeShared;
        qualityRight = qualityShared;
    }
public:
    // enable 5' cutting by quality
    bool enabledFront;
    // enable 3' cutting by quality
    bool enabledTail;
    // enable agressive cutting mode
    bool enabledRight;
    // the sliding window size
    int windowSizeShared;
    // the mean quality requirement
    int qualityShared;
    // the sliding window size for cutting by quality in 5'
    int windowSizeFront;
    // the mean quality requirement for cutting by quality in 5'
    int qualityFront;
    // the sliding window size for cutting by quality in 3'
    int windowSizeTail;
    // the mean quality requirement for cutting by quality in 3'
    int qualityTail;
    // the sliding window size for cutting by quality in aggressive mode
    int windowSizeRight;
    // the mean quality requirement for cutting by quality in aggressive mode
    int qualityRight;
};

class AdapterOptions {
public:
    AdapterOptions() {
        enabled = true;
        hasSeqR1 = false;
        hasSeqR2 = false;
        detectAdapterForPE = false;
        hasFasta = false;
    }

public:
    bool enabled;
    string sequence;
    string sequenceR2;
    string detectedAdapter1;
    string detectedAdapter2;
    vector<string> seqsInFasta;
    string fastaFile;
    bool hasSeqR1;
    bool hasSeqR2;
    bool hasFasta;
    bool detectAdapterForPE;
};

class TrimmingOptions {
public:
    TrimmingOptions() {
        front1 = 0;
        tail1 = 0;
        front2 = 0;
        tail2 = 0;
        maxLen1 = 0;
        maxLen2 = 0;
    }
public:
    // trimming first cycles for read1
    int front1;
    // trimming last cycles for read1
    int tail1;
    // trimming first cycles for read2
    int front2;
    // trimming last cycles for read2
    int tail2;
    // max length of read1
    int maxLen1;
    // max length of read2
    int maxLen2;
};

class QualityFilteringOptions {
public:
    QualityFilteringOptions() {
        enabled = true;
        // '0' = Q15
        qualifiedQual = '0';
        unqualifiedPercentLimit = 40;
        nBaseLimit = 5;
    }
public:
    // quality filter enabled
    bool enabled;
    // if a base's quality phred score < qualifiedPhred, then it's considered as a low_qual_base
    char qualifiedQual;
    // if low_qual_base_num > lowQualLimit, then discard this read
    int unqualifiedPercentLimit;
    // if n_base_number > nBaseLimit, then discard this read
    int nBaseLimit;
    // if average qual score < avgQualReq, then discard this read
    int avgQualReq;
};

class ReadLengthFilteringOptions {
    public:
        ReadLengthFilteringOptions() {
            enabled = false;
            requiredLength = 50;
            maxLength = 0;
        }
    public:
        // length filter enabled
        bool enabled;
        // if read_length < requiredLength, then this read is discard
        int requiredLength;
        // length limit, 0 for no limitation
        int maxLength;
};
enum Mode {
    MEM,
    GREEDY
};

enum CodonTable {
    codontable1,
    codontable2,
    codontable3,
    codontable4,
    codontable5,
    codontable6,
    codontable9,
    codontable10,
    codontable12,
    codontable13,
    codontable14,
    codontable16,
    codontable21,
    codontable22,
    codontable24,
    codontable26,
    codontable27,
    codontable29,
    codontable30,
    codontable31,
    codontable33
};
struct CommonSearchOptions{
    string fmi;
    bool SEG;
    bool useEvalue;
    double minEvalue;
    float lenper;
    unsigned int minFragLength;
    unsigned int misMatches;
    unsigned int minScore;
    unsigned int seedLength;
    size_t max_matches_SI;
    size_t max_match_ids;
    Mode mode;
    char db;
};
class TransSearchOptions {
public:
    TransSearchOptions() {
        codonTable = codontable1;
        comOptions.mode = GREEDY;
        comOptions.SEG = true;
        comOptions.useEvalue = false;
        comOptions.minEvalue = 0.01;
        comOptions.lenper = 0.9;
        comOptions.minFragLength = 11;
        comOptions.misMatches = 3;
        comOptions.minScore = 65;
        comOptions.seedLength = 7;
        comOptions.max_matches_SI = 100000;
        comOptions.max_match_ids = 100000;
        comOptions.db = 'p';
    }
    void reset2Default() {
    }

public:
    CommonSearchOptions comOptions;
    string tCodonTable;
    CodonTable codonTable;
};
class DNASearchOptions {
public:
    DNASearchOptions() {
        comOptions.mode = GREEDY;
        comOptions.SEG = true;
        comOptions.useEvalue = false;
        comOptions.minEvalue = 0.01;
        comOptions.lenper = 0.9;
        comOptions.minFragLength = 60;
        comOptions.misMatches = 6;
        comOptions.minScore = 50;
        comOptions.seedLength = 7;
        comOptions.max_matches_SI = 100000;
        comOptions.max_match_ids = 100000;
        comOptions.db = 'd';
    }

    void reset2Default() {
    }

public:
    CommonSearchOptions comOptions;
};
class HostSearchOptions {
public:
    HostSearchOptions() {
        comOptions.mode = GREEDY;
        comOptions.SEG = true;
        comOptions.useEvalue = false;
        comOptions.minEvalue = 0.01;
        comOptions.lenper = 0.9;
        comOptions.minFragLength = 90;
        comOptions.misMatches = 6;
        comOptions.minScore = 80;
        comOptions.seedLength = 7;
        comOptions.max_matches_SI = 1000;
        comOptions.max_match_ids = 1000;
        comOptions.db = 'h';
    }

    void reset2Default() {
    }

public:
    CommonSearchOptions comOptions;
};
class MarkerSearchOptions {
public:
    MarkerSearchOptions() {
        comOptions.mode = GREEDY;
        comOptions.SEG = true;
        comOptions.useEvalue = false;
        comOptions.minEvalue = 0.01;
        comOptions.lenper = 0.9;
        comOptions.minFragLength = 90;
        comOptions.misMatches = 6;
        comOptions.minScore = 80;
        comOptions.seedLength = 7;
        comOptions.max_matches_SI = 1000000;
        comOptions.max_match_ids = 1000000;
        comOptions.db = 'm';
    }

    void reset2Default() {
    }

public:
    CommonSearchOptions comOptions;
};
class HomoSearchOptions {
public:
    HomoSearchOptions() {
        mappedReads = 0;
        dnaReads = 0;
        proReads = 0;
        hostReads = 0;
        markerReads = 0;
        totReads = 0;
        totCleanReads = 0;
    }
    void reset2Default() {
        mappedReads = 0;
        dnaReads = 0;
        proReads = 0;
        hostReads = 0;
        markerReads = 0;
        totReads = 0;
        totCleanReads = 0;
    }

    double calculate(char readType = 't'){
        double mappedReadPer = 0.0;
        if(readType == 'd') {
            mappedReadPer = std::round((static_cast<double>(dnaReads * 100.0)/totCleanReads) * 10000.0) / 10000.0;
        } else if(readType == 'p'){
            mappedReadPer = std::round((static_cast<double>(proReads * 100.0)/totCleanReads) * 10000.0) / 10000.0;
        } else if(readType == 'h'){
            mappedReadPer = std::round((static_cast<double>(hostReads * 100.0)/totCleanReads) * 10000.0) / 10000.0;
        } else if(readType == 'm'){
            mappedReadPer = std::round((static_cast<double>(markerReads * 100.0)/totCleanReads) * 10000.0) / 10000.0;
        } else {
            mappedReadPer = std::round((static_cast<double>(mappedReads * 100.0)/totCleanReads) * 10000.0) / 10000.0;
        }
        return mappedReadPer;
    }

public:
    uint32 mappedReads;
    uint32 dnaReads;
    uint32 proReads;
    uint32 hostReads;
    uint32 markerReads;
    long totReads;
    long totCleanReads;
};

class Sample{
public:
    Sample(){
        prefix = "";
        ff = "";
        rf = "";
    }

    Sample(const std::string & p, const std::string & f, const std::string & r){
        prefix = p;
        ff = f;
        rf = r;
    };

public:
    std::string prefix;
    std::string ff;
    std::string rf;
};

struct Evaluation{
    bool supportEvaluation = false;
    bool disable_trim_poly_g = false;
};

class Options{
public:
    Options();
    ~Options();
    void init();
    bool isPaired();
    bool validate();
    bool adapterCuttingEnabled();
    bool polyXTrimmingEnabled();
    string getAdapter1();
    string getAdapter2();
    bool shallDetectAdapter(bool isR2 = false);
    void loadFastaAdapters();
    void deterCodonTable();
    void parseSampleTable();

public:
    // file name of read1 input
    string in1;
    // file name of read2 input
    string in2;
    // file name of read1 output
    string out1;
    // file name of read2 output
    string out2;
    string outdir;
    // output failed reads
    bool outFR;
    //file name of failed reads;
    string outFRFile;
    // json file
    string jsonFile;
    // html file
    string htmlFile;
    // html report title
    string reportTitle;
    // compression level
    int compression;
    // the input file is using phred64 quality scoring
    bool phred64;
    // do not rewrite existing files
    bool dontOverwrite;
    // read STDIN
    bool inputFromSTDIN;
    // write STDOUT
    bool outputToSTDOUT;
    // the input R1 file is interleaved
    bool interleavedInput;
    // only process first N reads
    int readsToProcess;
    // worker thread number
    int thread;
    // trimming options
    TrimmingOptions trim;
    // quality filtering options
    QualityFilteringOptions qualfilter;
    // length filtering options
    ReadLengthFilteringOptions lengthFilter;
    // adapter options
    AdapterOptions adapter;
    // options for quality cutting
    QualityCutOptions qualityCut;
    // options for base correction
    CorrectionOptions correction;
    // options for UMI
    UMIOptions umi;
    // 3' end polyG trimming, default for Illumina NextSeq/NovaSeq
    PolyGTrimmerOptions polyGTrim;
    // 3' end polyX trimming
    PolyXTrimmerOptions polyXTrim;
    int seqLen1;
    int seqLen2;
    // low complexity filtering
    LowComplexityFilterOptions complexityFilter;
    // options for duplication profiling
    DuplicationOptions duplicate;
    // options for duplication profiling
    int insertSizeMax;
    // overlap analysis threshold
    int overlapRequire;
    int overlapDiffLimit;
    int overlapDiffPercentLimit;
    // output debug information
    bool verbose;
    // the length of KMER, default is 25
    bool debug;

    std::string sampleTable;
    string prefix;
    std::vector<Sample> samVec;
    Evaluation mEvaluation;

    //merge overlapped PE read;
    bool mergerOverlappedPE;

    //FunTaxSeq file dir;
    string funtaxseqProgPath;
    //funtaxseq dir;
    string funtaxseqDir;

    string internalDBDir;
    //for internal database
    TransSearchOptions* mTransSearchOptions;
    HomoSearchOptions* mHomoSearchOptions;
    DNASearchOptions* mDNASearchOptions;
    HostSearchOptions* mHostSearchOptions;
    MarkerSearchOptions* mMarkerSearchOptions;
};

#endif
