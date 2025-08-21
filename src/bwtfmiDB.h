#ifndef BWTFMIDB_H
#define BWTFMIDB_H

#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <memory>
#include <thread>

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
    BwtFmiDB(Options * & opt);
    ~BwtFmiDB();

    void free_BWT();
    void free_FMI(FMI*& fmi);
    void free_suffixArray(suffixArray*& sa);
    void init(char db);

    BWT * bwt;
    FMI * fmi;
    AlphabetStruct * astruct;
    SegParameters * blast_seg_params;
    double db_length;
    bool transSearch;
    bool dnaSearch;
    bool hostSearch;
    bool markerSearch;
    std::string fmiFile;
    void print();

private:
    Options * mOptions;
    std::string database;
};

class BwtFmiDBPair{
public:
    BwtFmiDBPair(Options* & opt);
    void init();
    ~BwtFmiDBPair();
    bool transSearch;
    bool dnaSearch;
    bool hostSearch;
    bool markerSearch;

public:
    BwtFmiDB* tBwtfmiDB;
    BwtFmiDB* dBwtfmiDB;
    BwtFmiDB* hBwtfmiDB;
    BwtFmiDB* mBwtfmiDB;

private:
    Options * mOptions;
};
#endif /* BWTFMIDB_H */

