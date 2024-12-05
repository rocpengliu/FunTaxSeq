#include "homosearcher.h"

HomoSearcher::HomoSearcher(Options * & opt, BwtFmiDBPair* & bwtfmiDBPair) {
    mOptions = opt;
    match_ids.clear();
    locus = new std::string("");
    if (bwtfmiDBPair->transSearch){
        mTransSearcher = new TransSearcher(mOptions, bwtfmiDBPair->tBwtfmiDB);
        transSearch = true;
    } else {
        mTransSearcher = nullptr;
        transSearch = false;
    }

    if (bwtfmiDBPair->dnaSearch) {
        mDNASearcher = new DNASearcher(mOptions, bwtfmiDBPair->dBwtfmiDB);
        dnaSearch = true;
    } else {
        mDNASearcher = nullptr;
        dnaSearch = false;
    }
}

HomoSearcher::~HomoSearcher() {
    if (mTransSearcher) {
        delete mTransSearcher;
        mTransSearcher = nullptr;
    }

    if (mDNASearcher) {
        delete mDNASearcher;
        mDNASearcher = nullptr;
    }
    if(locus){
        delete locus;
        locus = nullptr;
    }
}

std::string* HomoSearcher::homoSearch(Read* & item) {
    locus->clear();
    match_ids.clear();
    if (mOptions->mDNASearchOptions->minFragLength < item->length())
        return locus;
    if (dnaSearch && transSearch) {
        match_ids = mDNASearcher->dnaSearchWuNeng(item);
        if (match_ids.empty()) {
            match_ids = mTransSearcher->transSearchWuKong(item);
        }
    } else if (dnaSearch) {
        match_ids = mDNASearcher->dnaSearchWuNeng(item);
    } else if (transSearch) {
        match_ids = mTransSearcher->transSearchWuKong(item);
    }
    if (!match_ids.empty()) {
       for(const auto &match : match_ids){
            if(match != nullptr){
                locus->append(match);
                locus->append(";");
            }
        }
    }
    return locus;
}

std::string* HomoSearcher::homoSearch(Read* & item1, Read* & item2) {
    locus->clear();
    match_ids.clear();
    if (mOptions->mDNASearchOptions->minFragLength <= item1->length() && mOptions->mDNASearchOptions->minFragLength <= item2->length()){
        if (dnaSearch && transSearch) {
            match_ids = mDNASearcher->dnaSearchWuNeng(item1, item2);
            if (match_ids.empty()) {
                match_ids = mTransSearcher->transSearchWuKong(item1, item2);
            }
        } else if (dnaSearch) {
            match_ids = mDNASearcher->dnaSearchWuNeng(item1, item2);
        } else if (transSearch) {
            match_ids = mTransSearcher->transSearchWuKong(item1, item2);
        }
    } else if (mOptions->mDNASearchOptions->minFragLength <= item1->length() && mOptions->mDNASearchOptions->minFragLength > item2->length()){
        if (dnaSearch && transSearch) {
            match_ids = mDNASearcher->dnaSearchWuNeng(item1);
            if (match_ids.empty()) {
                match_ids = mTransSearcher->transSearchWuKong(item1);
            }
        } else if (dnaSearch) {
            match_ids = mDNASearcher->dnaSearchWuNeng(item1);
        } else if (transSearch) {
            match_ids = mTransSearcher->transSearchWuKong(item1);
        }
    } else if (mOptions->mDNASearchOptions->minFragLength > item1->length() && mOptions->mDNASearchOptions->minFragLength <= item2->length()){
        if (dnaSearch && transSearch) {
            match_ids = mDNASearcher->dnaSearchWuNeng(item2);
            if (match_ids.empty()) {
                match_ids = mTransSearcher->transSearchWuKong(item2);
            }
        } else if (dnaSearch) {
            match_ids = mDNASearcher->dnaSearchWuNeng(item2);
        } else if (transSearch) {
            match_ids = mTransSearcher->transSearchWuKong(item2);
        }
    }
    if (!match_ids.empty()) {
        for(const auto &match : match_ids){
            if(match != nullptr){
                locus->append(match);
                locus->append(";");
            }
        }
    }
    return locus;
}
void HomoSearcher::clearMatchedIds(){
    for(auto it = match_ids.begin(); it != match_ids.end(); ++it) {
        if(*it) {
            delete *it;
        }
    }
    match_ids.clear();
}