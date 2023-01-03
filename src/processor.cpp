#include "processor.h"
#include "peprocessor.h"
#include "seprocessor.h"

Processor::Processor(Options* & opt){
    mOptions = opt;
    mBwtfmiDBPair = new BwtFmiDBPair(mOptions);
}

Processor::~Processor(){
    if(mBwtfmiDBPair){
        delete mBwtfmiDBPair;
        mBwtfmiDBPair = nullptr;
    }
}

bool Processor::process() {
    if(mOptions->isPaired()) {
        PairEndProcessor p(mOptions, mBwtfmiDBPair);
        p.process();
    } else {
        SingleEndProcessor p(mOptions, mBwtfmiDBPair);
        p.process();
    }

    return true;
}