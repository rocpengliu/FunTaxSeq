
#include "homosearcher.h"

HomoSearcher::HomoSearcher(Options * & opt, BwtFmiDBPair* & bwtfmiDBPair) {
    mOptions = opt;
    if(bwtfmiDBPair->transSearch){
        mTransSearcher = new TransSearcher(mOptions, bwtfmiDBPair->tBwtfmiDB);
        transSearch = true;
    } else {
        mTransSearcher = nullptr;
        transSearch = false;
    }
    
    if(bwtfmiDBPair->dnaSearch){
        mDNASearcher = new DNASearcher(mOptions, bwtfmiDBPair->dBwtfmiDB);
        dnaSearch = true;
    } else {
        mDNASearcher = nullptr;
        dnaSearch = false;
    }
}

HomoSearcher::~HomoSearcher() {
    if(mTransSearcher){
        delete mTransSearcher;
        mTransSearcher = nullptr;
    }
    
    if(mDNASearcher){
        delete mDNASearcher;
        mDNASearcher = nullptr;
    }
}

void HomoSearcher::homoSearch(Read* & item) {
    if (mOptions->mDNASearchOptions->minFragLength <= item->length()) {
        
        if(dnaSearch && transSearch){
            auto idd = mDNASearcher->dnaSearchWuNeng(item);

            if (idd == 0) {
                auto idp = mTransSearcher->transSearchWuKong(item);
                
                if(idp != 0 && mOptions->verbose){
                    cCout("protein", idp, 'b');
                }
                
            } else {
                
                if (mOptions->verbose) {
                     cCout("dna", idd, 'r');
                }
                
            }
            
        } else if(dnaSearch) {
            auto idd = mDNASearcher->dnaSearchWuNeng(item);
            if(idd != 0 && mOptions->verbose) {
                cCout("dna", idd, 'r');
                cCout(item->mSeq.mStr, 'r');
            }
        } else if(transSearch){
            auto idp = mTransSearcher->transSearchWuKong(item);
            if(idp != 0 && mOptions->verbose){
                 cCout("protein", idp, 'b');
            }
        }
    }
}

void HomoSearcher::homoSearch(Read* & item1, Read* & item2) {
    if (mOptions->mDNASearchOptions->minFragLength <= item1->length() && mOptions->mDNASearchOptions->minFragLength <= item2->length()) {
        if(dnaSearch && transSearch){
            auto idd = mDNASearcher->dnaSearchWuNeng(item1, item2);

            if (idd == 0) {
                auto idp = mTransSearcher->transSearchWuKong(item1, item2);

                if (idp != 0 && mOptions->verbose) {
                    cCout("protein", idp, 'b');
                }

            } else {

                if (mOptions->verbose) {
                    cCout("dna", idd, 'r');
                }

            }
            
            
        } else if(dnaSearch) {
            auto idd = mDNASearcher->dnaSearchWuNeng(item1, item2);
            if(idd != 0 && mOptions->verbose) {
                cCout("dna", idd, 'r');
            }
        } else if(transSearch){
            auto idp = mTransSearcher->transSearchWuKong(item1, item2);
            if(idp != 0 && mOptions->verbose){
                 cCout("protein", idp, 'b');
            }
        }

    } else if(mOptions->mDNASearchOptions->minFragLength <= item1->length() && mOptions->mDNASearchOptions->minFragLength > item2->length()){
        if(dnaSearch && transSearch){

            auto idd = mDNASearcher->dnaSearchWuNeng(item1);

            if (idd == 0) {
                auto idp = mTransSearcher->transSearchWuKong(item1);

                if (idp != 0 && mOptions->verbose) {
                    cCout("protein", idp, 'b');
                }

            } else {

                if (mOptions->verbose) {
                    cCout("dna", idd, 'r');
                }

            }
            
        } else if(dnaSearch) {
            auto idd = mDNASearcher->dnaSearchWuNeng(item1);
            if(idd != 0 && mOptions->verbose) {
                cCout("dna", idd, 'r');
            }
        } else if(transSearch){
            auto idp = mTransSearcher->transSearchWuKong(item1);
            if(idp != 0 && mOptions->verbose){
                 cCout("protein", idp, 'b');
            }
        }

    } else if(mOptions->mDNASearchOptions->minFragLength > item1->length() && mOptions->mDNASearchOptions->minFragLength <= item2->length()){
        if(dnaSearch && transSearch) {

            auto idd = mDNASearcher->dnaSearchWuNeng(item2);

            if (idd == 0) {
                auto idp = mTransSearcher->transSearchWuKong(item2);

                if (idp != 0 && mOptions->verbose) {
                    cCout("protein", idp, 'b');
                }

            } else {

                if (mOptions->verbose) {
                    cCout("dna", idd, 'r');
                }

            }
            
        } else if(dnaSearch) {
            auto idd = mDNASearcher->dnaSearchWuNeng(item2);
            if(idd != 0 && mOptions->verbose) {
                cCout("dna", idd, 'r');
            }
        } else if(transSearch){
            auto idp = mTransSearcher->transSearchWuKong(item2);
            if(idp != 0 && mOptions->verbose){
                 cCout("protein", idp, 'b');
            }
        }
    }
}

//void HomoSearcher::homoSearch(Read * item) {
//    if (mOptions->mDNASearchOptions->minFragLength <= item->length()) {
//        
//        if(dnaSearch && transSearch){
//            auto idd = mDNASearcher->dnaSearchWuNeng(item);
//            auto idp = mTransSearcher->transSearchWuKong(item);
//
//            if (mOptions->verbose) {
//                if (idd != 0 && idp != 0) {
//                    cCout("same", idd, 'g');
//                } else if (idd != 0 && idp == 0) {
//                    cCout("dna", idd, 'r');
//                } else if (idd == 0 && idp != 0) {
//                    cCout("protein", idp, 'b');
//                } else {
//
//                }
//            }
//            
//        } else if(dnaSearch) {
//            auto idd = mDNASearcher->dnaSearchWuNeng(item);
//            if(idd != 0 && mOptions->verbose) {
//                cCout("dna", idd, 'r');
//            }
//        } else if(transSearch){
//            auto idp = mTransSearcher->transSearchWuKong(item);
//            if(idp != 0 && mOptions->verbose){
//                 cCout("protein", idp, 'b');
//            }
//        }
//    }
//}
//
//void HomoSearcher::homoSearch(Read * item1, Read * item2) {
//    if (mOptions->mDNASearchOptions->minFragLength <= item1->length() && mOptions->mDNASearchOptions->minFragLength <= item2->length()) {
//        if(dnaSearch && transSearch){
//            auto idd = mDNASearcher->dnaSearchWuNeng(item1, item2);
//            auto idp = mTransSearcher->transSearchWuKong(item1, item2);
//
//            if (mOptions->verbose) {
//                if (idd != 0 && idp != 0) {
//                    cCout("same", idd, 'g');
//                } else if (idd != 0 && idp == 0) {
//                    cCout("dna", idd, 'r');
//                } else if (idd == 0 && idp != 0) {
//                    cCout("protein", idp, 'b');
//                } else {
//
//                }
//            }
//            
//        } else if(dnaSearch) {
//            auto idd = mDNASearcher->dnaSearchWuNeng(item1, item2);
//            if(idd != 0 && mOptions->verbose) {
//                cCout("dna", idd, 'r');
//            }
//        } else if(transSearch){
//            auto idp = mTransSearcher->transSearchWuKong(item1, item2);
//            if(idp != 0 && mOptions->verbose){
//                 cCout("protein", idp, 'b');
//            }
//        }
//
//    } else if(mOptions->mDNASearchOptions->minFragLength <= item1->length() && mOptions->mDNASearchOptions->minFragLength > item2->length()){
//        if(dnaSearch && transSearch){
//            auto idd = mDNASearcher->dnaSearchWuNeng(item1);
//            auto idp = mTransSearcher->transSearchWuKong(item1);
//
//            if (mOptions->verbose) {
//                if (idd != 0 && idp != 0) {
//                    cCout("same", idd, 'g');
//                } else if (idd != 0 && idp == 0) {
//                    cCout("dna", idd, 'r');
//                } else if (idd == 0 && idp != 0) {
//                    cCout("protein", idp, 'b');
//                } else {
//
//                }
//            }
//            
//        } else if(dnaSearch) {
//            auto idd = mDNASearcher->dnaSearchWuNeng(item1);
//            if(idd != 0 && mOptions->verbose) {
//                cCout("dna", idd, 'r');
//            }
//        } else if(transSearch){
//            auto idp = mTransSearcher->transSearchWuKong(item1);
//            if(idp != 0 && mOptions->verbose){
//                 cCout("protein", idp, 'b');
//            }
//        }
//
//    } else if(mOptions->mDNASearchOptions->minFragLength > item1->length() && mOptions->mDNASearchOptions->minFragLength <= item2->length()){
//        if(dnaSearch && transSearch){
//            auto idd = mDNASearcher->dnaSearchWuNeng(item2);
//            auto idp = mTransSearcher->transSearchWuKong(item2);
//
//            if (mOptions->verbose) {
//                if (idd != 0 && idp != 0) {
//                    cCout("same", idd, 'g');
//                } else if (idd != 0 && idp == 0) {
//                    cCout("dna", idd, 'r');
//                } else if (idd == 0 && idp != 0) {
//                    cCout("protein", idp, 'b');
//                } else {
//
//                }
//            }
//            
//        } else if(dnaSearch) {
//            auto idd = mDNASearcher->dnaSearchWuNeng(item2);
//            if(idd != 0 && mOptions->verbose) {
//                cCout("dna", idd, 'r');
//            }
//        } else if(transSearch){
//            auto idp = mTransSearcher->transSearchWuKong(item2);
//            if(idp != 0 && mOptions->verbose){
//                 cCout("protein", idp, 'b');
//            }
//        }
//    }
//}