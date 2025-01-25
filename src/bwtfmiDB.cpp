#include "bwtfmiDB.h"

BwtFmiDB::BwtFmiDB(Options * & opt) {
    mOptions = opt;
    fmiFile = "";
    transSearch = false;
    dnaSearch = false;
    hostSearch = false;
    markerSearch = false;
    astruct = NULL;
    blast_seg_params = NULL;
    database = "";
}

BwtFmiDB::~BwtFmiDB() {
    if (mOptions->verbose) {
        std::string msg = "start to free " + database + " database!";
        loginfo(msg);
    }
    free_BWT();
    free_AlphabetStruct(astruct);
    if (transSearch && blast_seg_params) {//mOptions->mTransSearchOptions->SEG
        //loginfo("Start to free protein blast_seg_params!");
        SegParametersFree(blast_seg_params);
        blast_seg_params = NULL;
    }
    if(dnaSearch && blast_seg_params) {//mOptions->mDNASearchOptions->SEG
        //loginfo("Start to free DNA blast_seg_params!");
        SegParametersFree(blast_seg_params);
        blast_seg_params = NULL;
    }
    if(hostSearch && blast_seg_params) {//mOptions->m
        //loginfo("Start to free host blast_seg_params!");
        SegParametersFree(blast_seg_params);
        blast_seg_params = NULL;
    }
    if(markerSearch && blast_seg_params) {//mOptions->mMarkerSearchOptions->SEG
        //loginfo("Start to free marker blast_seg_params!");
        SegParametersFree(blast_seg_params);
        blast_seg_params = NULL;
    }
    if (mOptions->verbose) {
        std::string msg = database + " database freed!";
        loginfo(msg);
    }
}

void BwtFmiDB::free_BWT() {
    //loginfo("Start to free bwt and fmi!");
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
        //loginfo("Start to free bwt!");
        free(bwt->bwt);
        bwt->bwt = NULL;
    }
    if (bwt->alphabet != NULL) {
        //loginfo("Start to free alphabet!");
        free(bwt->alphabet);
        bwt->alphabet = NULL;
    }
    // Finally, free the BWT structure itself
    free(bwt);
    bwt = NULL;
    //loginfo("Finished bwt and fmi freeing!");
}

void BwtFmiDB::free_FMI(FMI*& fmi) {
    //loginfo("Start to free fmi!");
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
    //loginfo("finish the fmi freeing!");
}

void BwtFmiDB::free_suffixArray(suffixArray*& sa) {
    //loginfo("Start to free suffixArray!");
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
    //loginfo("Finished suffixArray freeing!");
}

void BwtFmiDB::init(char db) {
    std::string msg = "";
    std::string dbfile = "";
    if (db == 'p' && !mOptions->mTransSearchOptions->comOptions.fmi.empty()) {
        if(file_exists(mOptions->mTransSearchOptions->comOptions.fmi)){
            dbfile = mOptions->mTransSearchOptions->comOptions.fmi;
            database = "protein";
        }
    } else if (db == 'd' && !mOptions->mDNASearchOptions->comOptions.fmi.empty()){
        if(file_exists(mOptions->mDNASearchOptions->comOptions.fmi)){
            dbfile = mOptions->mDNASearchOptions->comOptions.fmi;
            database = "DNA";
        }
    } else if(db == 'h' && !mOptions->mHostSearchOptions->comOptions.fmi.empty()){
        if(file_exists(mOptions->mHostSearchOptions->comOptions.fmi)){
            dbfile = mOptions->mHostSearchOptions->comOptions.fmi;
            database = "host";
        }
    } else if(db == 'm' && !mOptions->mMarkerSearchOptions->comOptions.fmi.empty()){
        if(file_exists(mOptions->mMarkerSearchOptions->comOptions.fmi)){
            dbfile = mOptions->mMarkerSearchOptions->comOptions.fmi;
            database = "marker";
        }
    } else {
        error_exit("must specify a database!");
    }
    msg.clear();
    msg = "Reading " + database + " BWT FMI index from file " + dbfile;
    if (mOptions->verbose) {
        loginfo(msg);
    }
    FILE * file = fopen(dbfile.c_str(), "r");
    bwt = readIndexes(file, db);
    fclose(file);
    fmi = bwt->f;
    db_length = (double) (bwt->len - bwt->nseq);
    astruct = alloc_AlphabetStruct(bwt->alphabet, 0, 0);
        //need to be conformed.
    if (db == 'p') {
        transSearch = true;
        if (mOptions->mTransSearchOptions->comOptions.SEG) {
            blast_seg_params = SegParametersNewAa(); //need to be conformed;
            blast_seg_params->overlaps = TRUE;
        }
    } else if (db == 'd'){
        dnaSearch = true;
        if (mOptions->mDNASearchOptions->comOptions.SEG) {
            blast_seg_params = SegParametersNewAa(); //need to be conformed;
            blast_seg_params->overlaps = TRUE;
        }
    } else if(db == 'h'){
        hostSearch = true;
        if (mOptions->mHostSearchOptions->comOptions.SEG) {
            blast_seg_params = SegParametersNewAa(); //need to be conformed;
            blast_seg_params->overlaps = TRUE;
        }
    } else if(db == 'm'){
        markerSearch = true;
        if (mOptions->mMarkerSearchOptions->comOptions.SEG) {
            blast_seg_params = SegParametersNewAa(); //need to be conformed;
            blast_seg_params->overlaps = TRUE;
        }
    } else {
        error_exit("must specify a database!");
    }
    if (mOptions->verbose) {
        msg = "finish BWT FMI initiation for " + database;
        loginfo(msg);
        msg.clear();
    };
}

void BwtFmiDB::print(){
    std::stringstream ss;
    std::string msg = "";
    ss << database << " BWT length: " << bwt->len << ", num seq: " << bwt->nseq << ", alphabet: " << bwt->alphabet;
    // "\nbwt->s->check: " << bwt->s->check <<
    // "\nbwt->s->chpt_exp: " << bwt->s->chpt_exp <<
    // "\nbwt->s->hash_step: " << bwt->s->hash_step <<
    // "\nbwt->s->mask: " << bwt->s->mask <<
    // "\nbwt->s->maxlength: " << bwt->s->maxlength <<
    // "\nbwt->s->nbytes: " << bwt->s->nbytes <<
    // "\nbwt->s->ncheck: " << bwt->s->ncheck <<
    // "\nbwt->s->pbits: " << bwt->s->pbits <<
    // "\nbwt->s->sbits: " << bwt->s->sbits <<
    // "\nfmi->N1: " << fmi->N1 <<
    // "\nfmi->N2: " << fmi->N2 <<
    // "\nfmi->alen: " << fmi->alen <<
    // "\nfmi->bwtlen: " << fmi->bwtlen <<
    // "db_length: " << db_length <<
    // "\n";
    msg = ss.str();
    loginfo(msg);
    msg.clear();
}

BwtFmiDBPair::BwtFmiDBPair(Options* & opt) {
    mOptions = opt;
    transSearch = false;
    dnaSearch = false;
    hostSearch = false;
    markerSearch = false;
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
    if(hBwtfmiDB) {
        delete hBwtfmiDB;
        hBwtfmiDB = nullptr;
    }
    if(mBwtfmiDB) {
        delete mBwtfmiDB;
        mBwtfmiDB = nullptr;
    }
}

void BwtFmiDBPair::init() {
    std::thread tThread([this](){
        if (mOptions->mTransSearchOptions->comOptions.fmi.empty()) {
            tBwtfmiDB = nullptr;
            transSearch = false;
        } else {
            tBwtfmiDB = new BwtFmiDB(mOptions);
            tBwtfmiDB->init('p');
            transSearch = true;
        }
    });

    std::thread dThread([this](){
        if (mOptions->mDNASearchOptions->comOptions.fmi.empty()) {
            dBwtfmiDB = nullptr;
            dnaSearch = false;
        } else {
            dBwtfmiDB = new BwtFmiDB(mOptions);
            dBwtfmiDB->init('d');
            dnaSearch = true;
        }
    });

    std::thread hThread([this](){
        if (mOptions->mHostSearchOptions->comOptions.fmi.empty()) {
            hBwtfmiDB = nullptr;
            hostSearch = false;
        } else {
            hBwtfmiDB = new BwtFmiDB(mOptions);
            hBwtfmiDB->init('h');
            hostSearch = true;
        }
    });

    std::thread mThread([this](){
        if(mOptions->mMarkerSearchOptions->comOptions.fmi.empty()){
            mBwtfmiDB = nullptr;
            markerSearch = false;
        } else {
            mBwtfmiDB = new BwtFmiDB(mOptions);
            mBwtfmiDB->init('m');
            markerSearch = true;
        }
    });

    if(tThread.joinable()){
        tThread.join();
    }
    if(dThread.joinable()){
        dThread.join();
    }
    if(hThread.joinable()){
        hThread.join();
    }

    if(mThread.joinable()){
        mThread.join();
    }

    if(hostSearch){
        hBwtfmiDB->print();
    }
    if(markerSearch){
        mBwtfmiDB->print();
    }
    if(dnaSearch){
        dBwtfmiDB->print();
    }
    if(transSearch){
        tBwtfmiDB->print();
    }
}