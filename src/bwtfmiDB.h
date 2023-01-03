#ifndef BWTFMIDB_H
#define BWTFMIDB_H

#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

#include "util.h"
#include "options.h"

#include "include/ncbi-blast+/algo/blast/core/blast_seg.h"
#include "include/ncbi-blast+/algo/blast/core/blast_filter.h"
#include "include/ncbi-blast+/algo/blast/core/blast_encoding.h"

extern "C" {
#include "bwt/fmi.h"
#include "bwt/bwt.h"
#include "bwt/sequence.h"
}
using namespace std;

class BwtFmiDB {
public:
    BwtFmiDB(Options * & opt, std::string db);
    ~BwtFmiDB();

    BWT * bwt;
    FMI * fmi;
    AlphabetStruct * astruct;
    SegParameters * blast_seg_params;
    double db_length;
    bool TransSearch;
    bool DNASearch;
    std::string fmiFile;
    
private:
    void init(std::string db);
    
private:
    Options * mOptions;
};

class BwtFmiDBPair{
public:
    BwtFmiDBPair(Options* & opt);
    ~BwtFmiDBPair();
    bool transSearch;
    bool dnaSearch;
public:
    BwtFmiDB* tBwtfmiDB;
    BwtFmiDB* dBwtfmiDB;
};
#endif /* BWTFMIDB_H */

