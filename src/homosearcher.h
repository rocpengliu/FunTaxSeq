#ifndef HOMOSEARCHER_H
#define HOMOSEARCHER_H

#include <valarray>

#include "options.h"
#include <memory>
#include "util.h"
#include "transsearcher.hpp"
#include "dnasearcher.h"
#include "read.h"

class HomoSearcher {
public:
    HomoSearcher(Options* & opt, BwtFmiDBPair* & bwtfmiDBPair);
    void homoSearch(Read* & item);
    void homoSearch(Read* & item1, Read* & item2);
    TransSearcher* getTransSearcher(){return mTransSearcher;};
    ~HomoSearcher();
    bool transSearch, dnaSearch;
private:
    Options * mOptions;
    TransSearcher * mTransSearcher;
    DNASearcher* mDNASearcher;
};

#endif /* HOMOSEARCHER_H */

