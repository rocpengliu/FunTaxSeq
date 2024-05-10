#include "processor.h"

Processor::Processor(Options* & opt){
    mOptions = opt;
    //mBwtfmiDBPair = bwtfmiDBPair;
}

Processor::~Processor(){
}

bool Processor::process(BwtFmiDBPair*& mBwtfmiDBPair) {
    if (mOptions->isPaired()){
        PairEndProcessor *p = new PairEndProcessor(mOptions, mBwtfmiDBPair);
        p->process();
        if (p){
            delete p;
            p = nullptr;
        }
    } else {
        SingleEndProcessor *p = new SingleEndProcessor(mOptions, mBwtfmiDBPair);
        p->process();
        if (p){
            delete p;
            p = nullptr;
        }
    }
    return true;
}