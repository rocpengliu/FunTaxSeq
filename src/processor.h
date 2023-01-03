#ifndef PROCESSOR_H
#define PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "options.h"
#include "bwtfmiDB.h"
#include <memory>

using namespace std;

class Processor{
public:
    Processor(Options* & opt);
    ~Processor();
    bool process();

private:
    Options* mOptions;
    BwtFmiDBPair* mBwtfmiDBPair;
};


#endif