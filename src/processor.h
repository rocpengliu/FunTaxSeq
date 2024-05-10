#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <memory>

#include "options.h"
#include "peprocessor.h"
#include "seprocessor.h"
#include "bwtfmiDB.h"

using namespace std;

class Processor{
public:
    Processor(Options* & opt);
    ~Processor();
    bool process(BwtFmiDBPair*& mBwtfmiDBPair);

private:
    Options* mOptions;
    BwtFmiDBPair* mBwtfmiDBPair;
    //PhyloTree* mPhyloTree;
};

#endif