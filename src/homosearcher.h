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
#include "util.h"
#include "transsearcher.hpp"
#include "dnasearcher.h"
#include "read.h"

class HomoSearcher {
public:
    HomoSearcher(Options* & opt, BwtFmiDBPair* & bwtfmiDBPair);
    std::string* homoSearch(Read* & item);
    std::string* homoSearch(Read* & item1, Read* & item2);
    ~HomoSearcher();
    bool transSearch, dnaSearch;
    HomoSearcher(const HomoSearcher & other) = delete;
    HomoSearcher & operator=(const HomoSearcher & other) = delete;
    std::string* locus;

private:
    void clearMatchedIds();
private:
    Options * mOptions;
    TransSearcher * mTransSearcher;
    DNASearcher* mDNASearcher;
    std::set<char*> match_ids;
};
#endif /* HOMOSEARCHER_H */