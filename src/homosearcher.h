#ifndef HOMOSEARCHER_H
#define HOMOSEARCHER_H

#include <valarray>
#include <memory>
#include <map>
#include <set>
#include <unordered_set>
#include <cstring>
#include <atomic>
#include <sstream>

#include "common.h"
#include "options.h"
#include "transsearcher.h"
#include "dnasearcher.h"
#include "read.h"

class HomoSearcher {
public:
    HomoSearcher(Options* & opt, BwtFmiDBPair* & bwtfmiDBPair);
    std::string* homoSearch(Read* & item, int & dnaReads, int & proReads, int & hostReads, int & markerReads);
    std::string* homoSearch(Read* & item1, Read* & item2, int & dnaReads, int & proReads, int & hostReads, int & markerReads);
    ~HomoSearcher();
    bool transSearch, dnaSearch, hostSearch, markerSearch;
    HomoSearcher(const HomoSearcher & other) = delete;
    HomoSearcher & operator=(const HomoSearcher & other) = delete;
    std::string* locus;

private:
    void clearMatchedIds();
    void postProcess();

private:
    Options * mOptions;
    TransSearcher* mTransSearcher;
    DNASearcher* mDNASearcher;
    DNASearcher* mHostSearcher;
    DNASearcher *mMarkerSearcher;
    std::set<char *> match_ids;
};
#endif /* HOMOSEARCHER_H */