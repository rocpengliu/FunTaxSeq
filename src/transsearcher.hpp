#ifndef TRANSSEARCHER_HPP
#define TRANSSEARCHER_HPP

#include <stdint.h>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <list>
#include <cmath>
#include <algorithm>
#include <mutex>
#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>
#include <string>
#include <cstring>
#include <climits>
#include <map>
#include <utility>
#include <functional>
#include <locale>
#include <stdio.h>
#include <cmath>
#include <valarray>
#include <memory>

#include "util.h"
#include "algo/blast/core/blast_seg.h"
#include "algo/blast/core/blast_filter.h"
#include "algo/blast/core/blast_encoding.h"
#include "read.h"
#include "options.h"
#include "fragment.h"
#include "options.h"
#include "bwtfmiDB.h"
#include "common.h"


extern "C" {
#include "bwt/bwt.h"
}

const double LN_2 = 0.6931471805;
const double LAMBDA = 0.3176;
const double LN_K = -2.009915479;

class TransSearcher {
protected:
    uint8_t codon_to_int(const char* codon);
    uint8_t revcomp_codon_to_int(const char* codon);

    uint8_t nuc2int[256];
    uint8_t compnuc2int[256];
    char codon2aa[256];
    uint8_t aa2int[256];

    std::map<char, std::vector<char>> blosum_subst;
    int8_t blosum62diag[20];
    int8_t b62[20][20];

    std::string translations[6];
    std::multimap<unsigned int, Fragment *, std::greater<unsigned int>> fragments;
    std::vector<SI *> best_matches_SI;
    std::vector<SI *> longest_matches_SI;
    std::vector<std::string> best_matches;
    std::vector<std::string> longest_fragments;
    
    unsigned int best_match_score = 0;
    double query_len = 0;
    uint32_t read_count = 0;
    uint32 uniq_mapped_reads = 0;
    uint32 multi_mapped_reads = 0;
    
    std::string readName;
    std::stringstream ss;

    void clearFragments();
    void clearMatchedIds();
    unsigned int calcScore(const std::string &);
    unsigned int calcScore(const std::string &, int);
    unsigned int calcScore(const std::string &, size_t, size_t, int);
    void addAllMismatchVariantsAtPosSI(const Fragment *, unsigned int, size_t, SI *); // used in Greedy mode
    Fragment * getNextFragment(unsigned int);
    void eval_match_scores(SI *si, Fragment *);
    void getAllFragmentsBits(const std::string& line);
    void getLongestFragmentsBits(const std::string& line);

protected:
    void classify_length();
    void classify_greedyblosum();
    void ids_from_SI(SI *);
    void ids_from_SI_recursive(SI *);
    void postProcess();
    std::set<char *> match_ids;

    Options * mOptions;   
    BwtFmiDB* tbwtfmiDB;
    
public:
    TransSearcher(Options * & opt, BwtFmiDB* & mBwtfmiDB);
    ~TransSearcher();
    std::set<char *>& transSearchWuKong(Read* & item);
    std::set<char *>& transSearchWuKong(Read* & item1, Read* & item2);
};

#endif /* TRANSSEARCHER_HPP */
