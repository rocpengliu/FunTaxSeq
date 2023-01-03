#ifndef SE_PROCESSOR_H
#define SE_PROCESSOR_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string>
#include <fstream>
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include <cstdlib>
#include <condition_variable>
#include <mutex>
#include <sstream>
#include <unordered_map>
#include <map>
#include <future>
#include <deque>
#include <time.h>
#include <iomanip>

#include "options.h"
#include "threadconfig.h"
#include "filter.h"
#include "umiprocessor.h"
#include "writerthread.h"
#include "duplicate.h"
#include "fastqreader.h"
#include "util.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "adaptertrimmer.h"
#include "polyx.h"
#include "read.h"
#include "common.h"
#include "bwtfmiDB.h"
#include <memory>


using namespace std;

struct ReadPack {
    Read** data;
    int count;
};

typedef struct ReadPack ReadPack;

struct ReadRepository {
    ReadPack** packBuffer;
    atomic_long readPos;
    atomic_long writePos;
};

typedef struct ReadRepository ReadRepository;

class SingleEndProcessor{
public:
    SingleEndProcessor(Options* & opt, BwtFmiDBPair* & bwtfmiDBPair);
    ~SingleEndProcessor();
    bool process();

private:
    bool processSingleEnd(ReadPack* pack, ThreadConfig* config);
    void initPackRepository();
    void destroyPackRepository();
    void producePack(ReadPack* pack);
    void consumePack(ThreadConfig* config);
    void producerTask();
    void consumerTask(ThreadConfig* config);
    void initConfig(ThreadConfig* config);
    void initOutput();
    void closeOutput();
    void writeTask(WriterThread* config);
private:
    Options* mOptions;
    ReadRepository mRepo;
    atomic_bool mProduceFinished;
    atomic_int mFinishedThreads;
    std::mutex mInputMtx;
    std::mutex mOutputMtx;
    std::mutex logMtx;
    Filter* mFilter;
    gzFile mZipFile;
    ofstream* mOutStream;
    UmiProcessor* mUmiProcessor;
    WriterThread* mLeftWriter;
    WriterThread* mFailedWriter;
    Duplicate* mDuplicate;
    BwtFmiDBPair* mBwtfmiDBPair;
};
#endif
