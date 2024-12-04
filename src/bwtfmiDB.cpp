#include "bwtfmiDB.h"

BwtFmiDB::BwtFmiDB(Options * & opt) {
    mOptions = opt;
    fmiFile = "";
    TransSearch = false;
    DNASearch = false;
    astruct = NULL;
    blast_seg_params = NULL;
}

BwtFmiDB::~BwtFmiDB() {
    free_BWT();
    free_AlphabetStruct(astruct);
    if (TransSearch && blast_seg_params) {//mOptions->mTransSearchOptions->SEG
        loginfo("Start to free protein blast_seg_params!");
        SegParametersFree(blast_seg_params);
        blast_seg_params = NULL;
    }
    if(DNASearch && blast_seg_params) {//mOptions->mDNASearchOptions->SEG
        loginfo("Start to free DNA blast_seg_params!");
        SegParametersFree(blast_seg_params);
        blast_seg_params = NULL;
    }
    if (mOptions->verbose) {
        if (TransSearch) {
            loginfo("Protein database freed!");
        }
        if (DNASearch) {
            loginfo("DNA database freed!");
        }
    }
}

void BwtFmiDB::free_BWT() {
    loginfo("Start to free bwt and fmi!");
    if (bwt == NULL)
        return; // Check if bwt is NULL
    if (bwt->f != NULL) {
        free_FMI(bwt->f);
        bwt->f = NULL;
    }
    if (bwt->s != NULL) {
        free_suffixArray(bwt->s);
        bwt->s = NULL;
    }
    // Free dynamically allocated members
    if (bwt->bwt != NULL) {
         loginfo("Start to free bwt!");
        free(bwt->bwt);
        bwt->bwt = NULL;
    }
    if (bwt->alphabet != NULL) {
         loginfo("Start to free alphabet!");
        free(bwt->alphabet);
        bwt->alphabet = NULL;
    }
    // Finally, free the BWT structure itself
    free(bwt);
    bwt = NULL;
    loginfo("Finished bwt and fmi freeing!");
}

void BwtFmiDB::free_FMI(FMI*& fmi) {
    loginfo("Start to free fmi!");
    if (fmi == NULL) return; // Check if fmi is NULL
    if (fmi->index1 != NULL) {
        for (int i = 0; i < fmi->N1; i++) {
            if (fmi->index1[i] != NULL) {
                free(fmi->index1[i]);
                fmi->index1[i] = NULL;
            }
        }
        free(fmi->index1);
        fmi->index1 = NULL;
    }
    if (fmi->index2 != NULL) {
        for (int i = 0; i < fmi->N2; i++) {
            if (fmi->index2[i] != NULL) {
                free(fmi->index2[i]);
                fmi->index2[i] = NULL;
            }
        }
        free(fmi->index2);
        fmi->index2 = NULL;
    }
    if (fmi->startLcode != NULL) {
        free(fmi->startLcode);
        fmi->startLcode = NULL;
    }
    
        // Free dynamically allocated members
    if (fmi->bwt != NULL) {
        free(fmi->bwt);
        fmi->bwt = NULL;
    }
    // Finally, free the FMI structure itself
    free(fmi);
    loginfo("finish the fmi freeing!");
}

void BwtFmiDB::free_suffixArray(suffixArray*& sa) {
    loginfo("Start to free suffixArray!");
    if (sa == NULL) return; // Check if sa is NULL
    // Free dynamically allocated members
    if (sa->sa != NULL) {
        free(sa->sa);
        sa->sa = NULL;
    }
    if (sa->seqTermOrder != NULL) {
        free(sa->seqTermOrder);
        sa->seqTermOrder = NULL;
    }
    if (sa->seqlengths != NULL) {
        free(sa->seqlengths);
        sa->seqlengths = NULL;
    }
    if (sa->hash != NULL) {
        for (int i = 0; i < sa->nseq; ++i) {
            SEQstruct *cur = sa->hash[i];
            recursive_free_SEQstruct(cur);
        }
        free(sa->hash);
        sa->hash = NULL;
    }
    if (sa->ids != NULL) {
        for (int i = 0; i < sa->nseq; ++i) {
            free(sa->ids[i]);
        }
        free(sa->ids);
        sa->ids = NULL;
    }
    //if (sa->seqstart != NULL) {
      //  free(sa->seqstart);
       // sa->seqstart = NULL;
   // }
    // Finally, free the suffixArray structure itself
    free(sa);
    sa = NULL;
    loginfo("Finished suffixArray freeing!");
}

void BwtFmiDB::init(std::string db) {
    std::stringstream ss;
    std::string msg = "";
    if (db == "protein" && !mOptions->mTransSearchOptions->tfmi.empty()) {
        ss.str("");
        msg.clear();
        msg = "Reading protein (trans search) BWT FMI index from file " + mOptions->mTransSearchOptions->tfmi;
        if (mOptions->verbose) {
            loginfo(msg);
        }
        FILE * file = fopen(mOptions->mTransSearchOptions->tfmi.c_str(), "r");
        bwt = readIndexes(file);
        fclose(file);
        fmi = bwt->f;
        if (mOptions->verbose) {
            ss << "Protein (trans search) BWT of length " << bwt->len << " has been read with " << bwt->nseq << " sequences, alphabet = " << bwt->alphabet <<
                    "\nbwt->s->check: " << bwt->s->check <<
                    "\nbwt->s->chpt_exp: " << bwt->s->chpt_exp <<
                    "\nbwt->s->hash_step: " << bwt->s->hash_step <<
                    "\nbwt->s->mask: " << bwt->s->mask <<
                    "\nbwt->s->maxlength: " << bwt->s->maxlength <<
                    "\nbwt->s->nbytes: " << bwt->s->nbytes <<
                    "\nbwt->s->nbytes: " << bwt->s->nbytes <<
                    "\nbwt->s->ncheck: " << bwt->s->ncheck <<
                    "\nbwt->s->pbits: " << bwt->s->pbits <<
                    "\nbwt->s->sbits: " << bwt->s->sbits <<
                    "\nfmi->N1: " << fmi->N1 <<
                    "\nfmi->N2: " << fmi->N2 <<
                    "\nfmi->alen: " << fmi->alen <<
                    "\nfmi->bwtlen: " << fmi->bwtlen;
            msg = ss.str();
            loginfo(msg);
        }
        db_length = (double) (bwt->len - bwt->nseq);
        msg.clear();
        if (mOptions->verbose) {
            msg = "Protein (trans search) double length is " + std::to_string(db_length);
            loginfo(msg);
        }
        astruct = alloc_AlphabetStruct(bwt->alphabet, 0, 0);
        //need to be conformed.
        if (mOptions->mTransSearchOptions->SEG) {
            blast_seg_params = SegParametersNewAa(); //need to be conformed;
            blast_seg_params->overlaps = TRUE;
        }
        TransSearch = true;
        if (mOptions->verbose) {
            loginfo("finish protein BWTFMI initiation");
        }
       // std::cin.get();
    }
    
    if (db == "DNA" && !mOptions->mDNASearchOptions->dfmi.empty()) {
        ss.str("");
        msg.clear();
        msg = "Reading DNA (DNA search) BWT FMI index from file " + mOptions->mDNASearchOptions->dfmi;
        if (mOptions->verbose) {
            loginfo(msg);
        }
        FILE * file = fopen(mOptions->mDNASearchOptions->dfmi.c_str(), "r");
        bwt = readIndexes(file);
        fclose(file);
        fmi = bwt->f;
        if (mOptions->verbose) {
            ss << "DNA (DNA) BWT of length " << bwt->len << " has been read with " << bwt->nseq << " sequences, alphabet = " << bwt->alphabet <<
                    "\nbwt->s->check: " << bwt->s->check <<
                    "\nbwt->s->chpt_exp: " << bwt->s->chpt_exp <<
                    "\nbwt->s->hash_step: " << bwt->s->hash_step <<
                    "\nbwt->s->mask: " << bwt->s->mask <<
                    "\nbwt->s->maxlength: " << bwt->s->maxlength <<
                    "\nbwt->s->nbytes: " << bwt->s->nbytes <<
                    "\nbwt->s->nbytes: " << bwt->s->nbytes <<
                    "\nbwt->s->ncheck: " << bwt->s->ncheck <<
                    "\nbwt->s->pbits: " << bwt->s->pbits <<
                    "\nbwt->s->sbits: " << bwt->s->sbits <<
                    "\nfmi->N1: " << fmi->N1 <<
                    "\nfmi->N2: " << fmi->N2 <<
                    "\nfmi->alen: " << fmi->alen <<
                    "\nfmi->bwtlen: " << fmi->bwtlen;
            msg = ss.str();
            loginfo(msg);
        }
        db_length = (double) (bwt->len - bwt->nseq);
        msg.clear();
        if (mOptions->verbose){
            msg = "DNA (DNA search) double length is " + std::to_string(db_length);
            loginfo(msg);
        }
        astruct = alloc_AlphabetStruct(bwt->alphabet, 0, 0);
        //need to be conformed.
        if (mOptions->mDNASearchOptions->SEG) {
            blast_seg_params = SegParametersNewAa(); //need to be conformed;
            blast_seg_params->overlaps = TRUE;
        }
        DNASearch = true;
        if (mOptions->verbose) {
            loginfo("finish DNA BWT FMI initiation");
        }
        //std::cin.get();
    }
}

BwtFmiDBPair::BwtFmiDBPair(Options* & opt) {
    mOptions = opt;
    init();
}

BwtFmiDBPair::~BwtFmiDBPair() {
    if (tBwtfmiDB) {
        delete tBwtfmiDB;
        tBwtfmiDB = nullptr;
    }
    if (dBwtfmiDB) {
        delete dBwtfmiDB;
        dBwtfmiDB = nullptr;
    }
}

void BwtFmiDBPair::init() {
    std::thread tThread([this](){
        if (mOptions->mTransSearchOptions->tfmi.empty()) {
            tBwtfmiDB = nullptr;
            transSearch = false;
        } else {
            tBwtfmiDB = new BwtFmiDB(mOptions);
            tBwtfmiDB->init("protein");
            transSearch = true;
        }
    });

    std::thread dThread([this](){
        if (mOptions->mDNASearchOptions->dfmi.empty()) {
            dBwtfmiDB = nullptr;
            dnaSearch = false;
        } else {
            dBwtfmiDB = new BwtFmiDB(mOptions);
            dBwtfmiDB->init("DNA");
            dnaSearch = true;
        }
    });

    if(tThread.joinable()){
        tThread.join();
    }
    if(dThread.joinable()){
        dThread.join();
    }
}