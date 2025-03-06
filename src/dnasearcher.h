#ifndef DNASEARCHER_H
#define DNASEARCHER_H

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

#include "util.h"
#include "algo/blast/core/blast_seg.h"
#include "algo/blast/core/blast_filter.h"
#include "algo/blast/core/blast_encoding.h"
#include "read.h"
#include "options.h"
#include "fragment.h"
#include "bwtfmiDB.h"
#include "common.h"

extern "C" {
#include "bwt/bwt.h"
}

// const double dLN_2 = 0.6931471805;
// const double dLAMBDA = 0.3176;
// const double dLN_K = -2.009915479;

class DNASearcher {
private:
    void classify_length();
    void classify_greedyblosum();
    void ids_from_SI(SI *);
    void ids_from_SI_recursive(SI *);
    void clearFragments();
    void clearMatchedIds();
    void getAllFragments(Read *&item);
    void postProcess();
    Fragment * getNextFragment(uint);
    void addAllMismatchVariantsAtPosSI(const Fragment *, uint, size_t, SI *); // used in Greedy mode
    uint calcScore(const std::string &);
    uint calcScore(const std::string &, int);
    uint calcScore(const std::string &, size_t, size_t, int);
    void eval_match_scores(SI *si, Fragment *);
private:
    uint8_t nuc2int[256];
    uint8_t compnuc2int[256];
    int8_t blosum62diag[4];
    int8_t b62[4][4];
    int query_len = 0;
    uint best_match_score = 0;
    std::multimap<uint, Fragment *, std::greater<uint>> fragments;
    std::vector<SI *> best_matches_SI;
    std::vector<SI *> longest_matches_SI;
    std::vector<std::string> best_matches;
    std::vector<std::string> longest_fragments;
    std::map<char, std::vector<char>> blosum_subst;
    //std::string readName;
    std::stringstream ss;
    uint min_match_length;

private:
    Options* mOptions;
    BwtFmiDB* mBwtfmiDB;
    CommonSearchOptions mCommonOptions;

public:
    DNASearcher(Options * & opt, BwtFmiDB* & bwtfmiDB, char db);
    ~DNASearcher();
    std::set<char *>& dnaSearchWuNeng(Read* & item);
    std::set<char *>& dnaSearchWuNeng(Read* & item1, Read* & item2);
    std::set<char *> match_ids;
};

#endif /* DNASEARCHER_H */

