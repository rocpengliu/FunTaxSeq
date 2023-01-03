#include "bwtfmiDB.h"

BwtFmiDB::BwtFmiDB(Options * & opt, std::string db) {
    mOptions = opt;
    fmiFile = "";
    TransSearch = false;
    DNASearch = false;
    init(db);
}

BwtFmiDB::~BwtFmiDB() {
    if (astruct->trans) free(astruct->trans);
    if (astruct->a) free(astruct->a);
    if (astruct) free(astruct);
    if (mOptions->mTransSearchOptions->SEG && TransSearch) {
        SegParametersFree(blast_seg_params);
    }
    if (mOptions->mDNASearchOptions->SEG && DNASearch) {
        SegParametersFree(blast_seg_params);
    }
    if (mOptions->verbose) {
        if(TransSearch){
            mOptions->longlog ? loginfolong("Protein database deleted!") : loginfo("Protein database deleted!");
        }

        if (DNASearch) {
            mOptions->longlog ? loginfolong("DNA database deleted!") : loginfo("DNA database deleted!");
        }
    }
}

void BwtFmiDB::init(std::string db) {

    std::string msg = "";
    fmiFile = "";
    if (db == "protein" && !mOptions->mTransSearchOptions->tfmi.empty()) {
        fmiFile = mOptions->mTransSearchOptions->tfmi;
        
        msg = "Reading protein (trans search) BWT FMI index from file " + mOptions->mTransSearchOptions->tfmi;
        if (mOptions->verbose) {
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
        }

        FILE * file = fopen(fmiFile.c_str(), "r");
        bwt = readIndexes(file);
        TransSearch = true;
        fclose(file);
        fmi = bwt->f;
        
        std::stringstream ss;
        std::string msg = "";
        if (mOptions->verbose) {
            ss << "Protein (trans search) BWT of length "  << bwt->len << " has been read with " << bwt->nseq << " sequences, alphabet = " << bwt->alphabet <<
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
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
        }

        db_length = (double) (bwt->len - bwt->nseq);
        
        msg = "";
        if (mOptions->verbose) {
            msg = "Protein (trans search) double length is " + std::to_string(db_length);
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
        }

        astruct = alloc_AlphabetStruct(bwt->alphabet, 0, 0);

        //need to be conformed.
        if (mOptions->mTransSearchOptions->SEG) {
            blast_seg_params = SegParametersNewAa(); //need to be conformed;
            blast_seg_params->overlaps = TRUE;
        }

        if (mOptions->verbose) {
            mOptions->longlog ? loginfolong("finish protein BWTFMI initiation") : loginfo("finish protein BWTFMI initiation");
        }
    
    } else if(db == "DNA" && !mOptions->mDNASearchOptions->dfmi.empty()) {
        fmiFile = mOptions->mDNASearchOptions->dfmi;
        msg = "Reading DNA (DNA search) BWT FMI index from file " + mOptions->mDNASearchOptions->dfmi;
        
        if (mOptions->verbose) {
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
        }

        FILE * file = fopen(fmiFile.c_str(), "r");
        bwt = readIndexes(file);
        DNASearch = true;
        fclose(file);
        fmi = bwt->f;
        
        std::stringstream ss;
        std::string msg = "";
        if (mOptions->verbose) {
            ss << "DNA (DNA) BWT of length "  << bwt->len << " has been read with " << bwt->nseq << " sequences, alphabet = " << bwt->alphabet <<
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
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
        }

        db_length = (double) (bwt->len - bwt->nseq);
        
        msg = "";
        if (mOptions->verbose) {
            msg = "DNA (DNA search) double length is " + std::to_string(db_length);
            mOptions->longlog ? loginfolong(msg) : loginfo(msg);
        }

        astruct = alloc_AlphabetStruct(bwt->alphabet, 0, 0);

        //need to be conformed.
        if (mOptions->mDNASearchOptions->SEG) {
            blast_seg_params = SegParametersNewAa(); //need to be conformed;
            blast_seg_params->overlaps = TRUE;
        }
        DNASearch = true;
        if (mOptions->verbose) {
            mOptions->longlog ? loginfolong("finish DNA BWT FMI initiation") : loginfo("finish DNA BWT FMI initiation");
        }
    } else {
        error_exit("invalide database, you have to supply the correct protein and/or DNA database!");
    }
    
}

BwtFmiDBPair::BwtFmiDBPair(Options* & opt) {

    if (opt->mTransSearchOptions->tfmi.empty()) {
        tBwtfmiDB = nullptr;
        transSearch = false;
    } else {
        tBwtfmiDB = new BwtFmiDB(opt, "protein");
        transSearch = true;
    }

    if (opt->mDNASearchOptions->dfmi.empty()) {
        dBwtfmiDB = nullptr;
        dnaSearch = false;
    } else {
        dBwtfmiDB = new BwtFmiDB(opt, "DNA");
        dnaSearch = true;
    }

}

BwtFmiDBPair::~BwtFmiDBPair(){
    if(tBwtfmiDB){
        delete tBwtfmiDB;
        tBwtfmiDB = nullptr;
    }
    if(dBwtfmiDB){
        delete dBwtfmiDB;
        dBwtfmiDB = nullptr;
    }
}