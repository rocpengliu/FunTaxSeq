#include "homosearcher.h"
#include "util.h"

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
        mDNASearcher = new DNASearcher(mOptions, bwtfmiDBPair->dBwtfmiDB, 'd');
        dnaSearch = true;
    } else {
        mDNASearcher = nullptr;
        dnaSearch = false;
    }
    if(bwtfmiDBPair->hostSearch){
        mHostSearcher = new DNASearcher(mOptions, bwtfmiDBPair->hBwtfmiDB, 'h');
        hostSearch = true;
    } else{
        mHostSearcher = nullptr;
        hostSearch = false;
    }
    if(bwtfmiDBPair->markerSearch) {
        mMarkerSearcher = new DNASearcher(mOptions, bwtfmiDBPair->mBwtfmiDB, 'm');
        markerSearch = true;
    } else {
        mMarkerSearcher = nullptr;
        markerSearch = false;
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
    if(mHostSearcher){
        delete mHostSearcher;
        mHostSearcher = nullptr;
    }
    if(mMarkerSearcher){
        delete mMarkerSearcher;
        mMarkerSearcher = nullptr;
    }
    if(locus){
        delete locus;
        locus = nullptr;
    }
}

void HomoSearcher::postProcess(char type){
    if (!match_ids.empty()) {
        for(const auto & match : match_ids){
            if(match != nullptr){
                //cCout("match id is", match, 'b');
                locus->append(match);
                locus->append(";");
                if(type == 'd'){

                } else if(type == 'p'){

                } else {

                }
            }
        }
    }
}

std::string* HomoSearcher::homoSearch(Read* & item, int & dnaReads, int & proReads, int & hostReads, int & markerReads) {
    locus->clear();
    match_ids.clear();
    if(hostSearch){
        if(mOptions->mHostSearchOptions->comOptions.minFragLength <= item->length()){
            match_ids = mHostSearcher->dnaSearchWuNeng(item);
            if(!match_ids.empty()){
                ++hostReads;
                locus->append(trimName2(item->mName) + "\thost\t");
                postProcess('h');
                return locus;
            }
        }
    }

    if(markerSearch){
        if(mOptions->mMarkerSearchOptions->comOptions.minFragLength <= item->length()){
            match_ids = mMarkerSearcher->dnaSearchWuNeng(item);
            if(!match_ids.empty()){
                ++markerReads;
                locus->append(trimName2(item->mName) + "\tmarker\t");
                postProcess('m');
                return locus;
            }
        }
    }

    if (dnaSearch) {
        if(mOptions->mDNASearchOptions->comOptions.minFragLength <= item->length()){
            match_ids = mDNASearcher->dnaSearchWuNeng(item);
            if (!match_ids.empty()){
                ++dnaReads;
                locus->append(trimName2(item->mName) + "\tdna\t");
                postProcess('d');
                return locus;
            }
        }
    }

    if (transSearch) {
        if(mOptions->mTransSearchOptions->comOptions.minFragLength * 3 <= item->length()){
            match_ids = mTransSearcher->transSearchWuKong(item);
            if(!match_ids.empty()){
                ++proReads;
                locus->append(trimName2(item->mName) + "\tpro\t");
                postProcess('m');
                return locus;
            }
        }
    }
    return locus;
}

std::string* HomoSearcher::homoSearch(Read* & item1, Read* & item2, int & dnaReads, int & proReads, int & hostReads, int & markerReads) {
    locus->clear();
    match_ids.clear();

    if(hostSearch){
        if(mOptions->mHostSearchOptions->comOptions.minFragLength <= min(item1->length(), item2->length())){
            match_ids = mHostSearcher->dnaSearchWuNeng(item1, item2);
        } else if(mOptions->mHostSearchOptions->comOptions.minFragLength <= item1->length()){
            match_ids = mHostSearcher->dnaSearchWuNeng(item1);
        } else if(mOptions->mHostSearchOptions->comOptions.minFragLength <= item2->length()){
            match_ids = mHostSearcher->dnaSearchWuNeng(item2);
        }

        if(!match_ids.empty()){
            ++hostReads;
            locus->append(trimName2(item1->mName) + "\thost\t");
            postProcess('h');
            return locus;
        }
    }

    if(markerSearch){
        if(mOptions->mMarkerSearchOptions->comOptions.minFragLength <= min(item1->length(), item2->length())){
            match_ids = mMarkerSearcher->dnaSearchWuNeng(item1, item2);
        } else if(mOptions->mMarkerSearchOptions->comOptions.minFragLength <= item1->length()){
            match_ids = mMarkerSearcher->dnaSearchWuNeng(item1);
        } else if(mOptions->mMarkerSearchOptions->comOptions.minFragLength <= item2->length()){
            match_ids = mMarkerSearcher->dnaSearchWuNeng(item2);
        }

        if(!match_ids.empty()){
            ++markerReads;
            locus->append(trimName2(item1->mName) + "\tmarker\t");
            postProcess('m');
            return locus;
        }
    }

    if (dnaSearch){
        if(mOptions->mDNASearchOptions->comOptions.minFragLength <= min(item1->length(), item2->length())){
            match_ids = mDNASearcher->dnaSearchWuNeng(item1, item2);
        } else if(mOptions->mDNASearchOptions->comOptions.minFragLength <= item1->length()){
            match_ids = mDNASearcher->dnaSearchWuNeng(item1);
        } else if(mOptions->mDNASearchOptions->comOptions.minFragLength <= item2->length()){
            match_ids = mDNASearcher->dnaSearchWuNeng(item2);
        }

        if(!match_ids.empty()){
            ++dnaReads;
            locus->append(trimName2(item1->mName) + "\tdna\t");
            postProcess('d');
            return locus;
        }
    }

    if(transSearch){
        if(mOptions->mTransSearchOptions->comOptions.minFragLength * 3 <= min(item1->length(), item2->length())){
            match_ids = mTransSearcher->transSearchWuKong(item1, item2);
        } else if(mOptions->mTransSearchOptions->comOptions.minFragLength * 3 <= item1->length()){
            match_ids = mTransSearcher->transSearchWuKong(item1);
        } else if(mOptions->mTransSearchOptions->comOptions.minFragLength * 3 <= item2->length()){
            match_ids = mTransSearcher->transSearchWuKong(item2);
        }

        if(!match_ids.empty()){
            ++proReads;
            locus->append(trimName2(item1->mName) + "\tpro\t");
            postProcess('p');
            return locus;
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