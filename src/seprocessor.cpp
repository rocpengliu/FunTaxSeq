#include "seprocessor.h"
#include "fastqreader.h"
#include <iostream>
#include <unistd.h>
#include <functional>
#include <thread>
#include <memory.h>
#include "util.h"
#include "jsonreporter.h"
#include "htmlreporter.h"
#include "adaptertrimmer.h"
#include "polyx.h"

SingleEndProcessor::SingleEndProcessor(Options* opt, BwtFmiDBPair* & bwtfmiDBPair){
    mOptions = opt;
    mProduceFinished = false;
    mFinishedThreads = 0;
    mFilter = new Filter(opt);
    mZipFile = NULL;
    mUmiProcessor = new UmiProcessor(opt);
    mLeftWriter =  NULL;
    mFailedWriter = NULL;

    mDuplicate = NULL;
    if(mOptions->duplicate.enabled) {
        mDuplicate = new Duplicate(mOptions);
    }
    mBwtfmiDBPair = bwtfmiDBPair;
    //mPhyloTree = phyloTree;
}

SingleEndProcessor::~SingleEndProcessor() {
    if(mFilter){
        delete mFilter;
        mFilter = NULL;
    }
    if(mDuplicate) {
        delete mDuplicate;
        mDuplicate = NULL;
    }
    
    if(mUmiProcessor){
        delete mUmiProcessor;
        mUmiProcessor = NULL;
    }
}

void SingleEndProcessor::initOutput() {
    if(!mOptions->out1.empty())
        mLeftWriter = new WriterThread(mOptions, mOptions->out1);
    if (!mOptions->outFRFile.empty()){
        mFailedWriter = new WriterThread(mOptions, mOptions->outFRFile);
    }
}

void SingleEndProcessor::closeOutput() {
    if(mLeftWriter) {
        delete mLeftWriter;
        mLeftWriter = NULL;
    }
    
    if (mFailedWriter) {
        delete mFailedWriter;
        mFailedWriter = NULL;
    }
}

void SingleEndProcessor::initConfig(ThreadConfig* config) {
    if(mOptions->out1.empty())
        return;
}

bool SingleEndProcessor::process(){
    initOutput();
    initPackRepository();
    std::thread producer(std::bind(&SingleEndProcessor::producerTask, this));

    //TODO: get the correct cycles
    int cycle = 151;
    ThreadConfig** configs = new ThreadConfig*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        configs[t] = new ThreadConfig(mOptions, mBwtfmiDBPair, t, false);
        initConfig(configs[t]);
    }

    std::thread** threads = new thread*[mOptions->thread];
    for(int t=0; t<mOptions->thread; t++){
        threads[t] = new std::thread(std::bind(&SingleEndProcessor::consumerTask, this, configs[t]));
    }

    std::thread* leftWriterThread = NULL;
    if(mLeftWriter)
        leftWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mLeftWriter));

    std::thread* failedWriterThread = NULL;
    if (mFailedWriter) {
        failedWriterThread = new std::thread(std::bind(&SingleEndProcessor::writeTask, this, mFailedWriter));
    }

    producer.join();
    for(int t=0; t<mOptions->thread; t++){
        threads[t]->join();
    }

    if (leftWriterThread)
        leftWriterThread->join();
    
    if (failedWriterThread)
        failedWriterThread->join();

    destroyPackRepository();
    
    if(mOptions->verbose)
        loginfo("start to generate reports\n");

    // merge stats and read filter results
    vector<Stats*> preStats;
    vector<Stats*> postStats;
    vector<FilterResult*> filterResults;

    for(int t=0; t<mOptions->thread; t++){
        preStats.push_back(configs[t]->getPreStats1());
        postStats.push_back(configs[t]->getPostStats1());
        filterResults.push_back(configs[t]->getFilterResult());
    }
    
    Stats* finalPreStats = Stats::merge(preStats);
    Stats* finalPostStats = Stats::merge(postStats);
    FilterResult* finalFilterResult = FilterResult::merge(filterResults);

    int* dupHist = NULL;
    double* dupMeanTlen = NULL;
    double* dupMeanGC = NULL;
    double dupRate = 0.0;
    if(mOptions->duplicate.enabled) {
        dupHist = new int[mOptions->duplicate.histSize];
        memset(dupHist, 0, sizeof(int) * mOptions->duplicate.histSize);
        dupMeanGC = new double[mOptions->duplicate.histSize];
        memset(dupMeanGC, 0, sizeof(double) * mOptions->duplicate.histSize);
        dupRate = mDuplicate->statAll(dupHist, dupMeanGC, mOptions->duplicate.histSize);
        cerr << endl;
        cerr << "Duplication rate (may be overestimated since this is SE data): " << dupRate * 100.0 << "%" << endl;
    }

    JsonReporter jr(mOptions);
    jr.setDupHist(dupHist, dupMeanGC, dupRate);
    jr.report(finalFilterResult, finalPreStats, finalPostStats);
    cerr << "Finished Json report" << endl;
    // make HTML report
    HtmlReporter hr(mOptions);
    hr.setDupHist(dupHist, dupMeanGC, dupRate);
    hr.report(finalFilterResult, finalPreStats, finalPostStats);
    cerr << "Finished Html report" << endl;

    // clean up
    for(int t=0; t<mOptions->thread; t++){
        delete threads[t];
        threads[t] = NULL;
        delete configs[t];
        configs[t] = NULL;
    }

    delete finalPreStats;
    delete finalPostStats;
    delete finalFilterResult;

    if(mOptions->duplicate.enabled) {
        delete[] dupHist;
        delete[] dupMeanGC;
    }

    delete[] threads;
    delete[] configs;

    if(leftWriterThread)
        delete leftWriterThread;

    if (failedWriterThread)
        delete failedWriterThread;

    closeOutput();

    return true;
}

bool SingleEndProcessor::processSingleEnd(ReadPack* pack, ThreadConfig* config){
    string outstr = "";
    string failedOutput = "";
    string locus = "";
    int readPassed = 0;
    int dnaReads = 0;
    int proReads = 0;
    int hostReads = 0;
    int markerReads = 0;
    for(int p=0;p<pack->count;p++){

        // original read1
        Read* or1 = pack->data[p];

        // stats the original read before trimming
        config->getPreStats1()->statRead(or1);

        // handling the duplication profiling
        if(mDuplicate)
            mDuplicate->statRead(or1);

        // umi processing
        if(mOptions->umi.enabled)
            mUmiProcessor->process(or1);

        int frontTrimmed = 0;
        // trim in head and tail, and apply quality cut in sliding window
        //Read* r1 = mFilter->trimAndCut(or1, mOptions->trim.front1, mOptions->trim.tail1, frontTrimmed);
        Read* r1 = mFilter->trimAndCut(or1, 0, 0, frontTrimmed);
        if(r1 != NULL) {
            if(mOptions->polyGTrim.enabled)
                PolyX::trimPolyG(r1, config->getFilterResult(), mOptions->polyGTrim.minLen);
        } 

        if(r1 != NULL && mOptions->adapter.enabled){
            bool trimmed = false;
            if(mOptions->adapter.hasSeqR1)
                trimmed = AdapterTrimmer::trimBySequence(r1, config->getFilterResult(), mOptions->adapter.sequence, false);
            bool incTrimmedCounter = !trimmed;
            if(mOptions->adapter.hasFasta) {
                AdapterTrimmer::trimByMultiSequences(r1, config->getFilterResult(), mOptions->adapter.seqsInFasta, false, incTrimmedCounter);
            }
        }

        if(r1 != NULL) {
            if(mOptions->polyXTrim.enabled)
                PolyX::trimPolyX(r1, config->getFilterResult(), mOptions->polyXTrim.minLen);
        }

        if(r1 != NULL) {
            if( mOptions->trim.maxLen1 > 0 && mOptions->trim.maxLen1 < r1->length())
                r1->resize(mOptions->trim.maxLen1);
        }

        int result = mFilter->passFilter(r1);

        config->addFilterResult(result, 1);

        if( r1 != NULL &&  result == PASS_FILTER) {
            locus.clear();
            locus = *(config->getHomoSearcher()->homoSearch(r1, dnaReads, proReads, hostReads, markerReads));
            if (!locus.empty()) {
                failedOutput += (locus + "\n");
            } else{
                outstr += r1->toString();
            }
            // stats the read after filtering
            config->getPostStats1()->statRead(r1);
            ++readPassed;
        }

        // if no trimming applied, r1 should be identical to or1
        if(r1 != or1 && r1 != NULL){
            delete r1;
            r1 = NULL;
        } else if(r1 == NULL){
            delete or1;
            or1 = NULL;
        }
    }
    // if splitting output, then no lock is need since different threads write different files
    mOutputMtx.lock();
    mOptions->mHomoSearchOptions->dnaReads += dnaReads;
    mOptions->mHomoSearchOptions->proReads += proReads;
    mOptions->mHomoSearchOptions->hostReads += hostReads;
    mOptions->mHomoSearchOptions->markerReads += markerReads;
    mOptions->mHomoSearchOptions->mappedReads += (dnaReads + proReads + hostReads + markerReads);
    if (mOptions->mHomoSearchOptions->mappedReads % 1000 == 0){
        cerr << "mapped " << (mOptions->mHomoSearchOptions->mappedReads / 1000) << "K reads("
            << (mOptions->mHomoSearchOptions->hostReads / 1000) << "K h|"
            << (mOptions->mHomoSearchOptions->markerReads / 1000) << "K m|"
            << (mOptions->mHomoSearchOptions->dnaReads / 1000) << "K d|"
            << (mOptions->mHomoSearchOptions->proReads / 1000) << "K p"
            << ")\n";
    }
    if(mOptions->outputToSTDOUT) {
        fwrite(outstr.c_str(), 1, outstr.length(), stdout);
    }

    if(mLeftWriter) {
        char* ldata = new char[outstr.size()];
        memcpy(ldata, outstr.c_str(), outstr.size());
        mLeftWriter->input(ldata, outstr.size());
    }

    if (mFailedWriter && !failedOutput.empty()) {
        char *fdata = new char[failedOutput.size()];
        memcpy(fdata, failedOutput.c_str(), failedOutput.size());
        mFailedWriter->input(fdata, failedOutput.size());
    }
    mOutputMtx.unlock();
    config->markProcessed(pack->count);
    delete[] pack->data;
    delete pack;
    return true;
}

void SingleEndProcessor::initPackRepository() {
    mRepo.packBuffer = new ReadPack*[PACK_NUM_LIMIT];
    memset(mRepo.packBuffer, 0, sizeof(ReadPack*)*PACK_NUM_LIMIT);
    mRepo.writePos = 0;
    mRepo.readPos = 0;
}

void SingleEndProcessor::destroyPackRepository() {
    delete[] mRepo.packBuffer;
    mRepo.packBuffer = NULL;
}

void SingleEndProcessor::producePack(ReadPack* pack){
    mRepo.packBuffer[mRepo.writePos] = pack;
    mRepo.writePos++;
}

void SingleEndProcessor::consumePack(ThreadConfig* config){
    ReadPack* data;
    mInputMtx.lock();
    while(mRepo.writePos <= mRepo.readPos) {
        usleep(1000);
        if(mProduceFinished) {
            mInputMtx.unlock();
            return;
        }
    }
    data = mRepo.packBuffer[mRepo.readPos];
    mRepo.readPos++;
    mInputMtx.unlock();
    processSingleEnd(data, config);

}

void SingleEndProcessor::producerTask(){
    if(mOptions->verbose)
        loginfo("start to load data");
    long lastReported = 0;
    int slept = 0;
    long readNum = 0;
    bool splitSizeReEvaluated = false;
    Read** data = new Read*[PACK_SIZE];
    memset(data, 0, sizeof(Read*)*PACK_SIZE);
    FastqReader reader(mOptions->in1, true, mOptions->phred64);
    int count=0;
    bool needToBreak = false;
    while(true){
        Read* read = reader.read();
        // TODO: put needToBreak here is just a WAR for resolve some unidentified dead lock issue 
        if(!read || needToBreak){
            // the last pack
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            data = NULL;
            if(read) {
                delete read;
                read = NULL;
            }
            break;
        }
        data[count] = read;
        ++count;
        // configured to process only first N reads
        if(mOptions->readsToProcess >0 && count + readNum >= mOptions->readsToProcess) {
            needToBreak = true;
        }
        if(mOptions->verbose && count + readNum >= lastReported + 1000000) {
            lastReported = count + readNum;
            string msg = "loaded " + to_string((lastReported/1000000)) + "M reads";
            loginfo(msg);
        }
        // a full pack
        if(count == PACK_SIZE || needToBreak){
            ReadPack* pack = new ReadPack;
            pack->data = data;
            pack->count = count;
            producePack(pack);
            //re-initialize data for next pack
            data = new Read*[PACK_SIZE];
            memset(data, 0, sizeof(Read*)*PACK_SIZE);
            // if the consumer/writer is far behind this producer/reader, sleep and wait to limit memory usage
            while(mRepo.writePos - mRepo.readPos > PACK_IN_MEM_LIMIT){
                ++slept;
                usleep(100);
            }
            readNum += count;
            // if the writer threads are far behind this producer, sleep and wait
            // check this only when necessary
            if(readNum % (PACK_SIZE * PACK_IN_MEM_LIMIT) == 0 && mLeftWriter) {
                while(mLeftWriter->bufferLength() > PACK_IN_MEM_LIMIT || (mFailedWriter && mFailedWriter->bufferLength() > PACK_IN_MEM_LIMIT)) {
                    ++slept;
                    usleep(1000);
                }
            }
            // reset count to 0
            count = 0;
        }
    }

    mProduceFinished = true;
    if(mOptions->verbose)
        loginfo("all reads loaded, start to monitor thread status");

    // if the last data initialized is not used, free it
     if(data != NULL){
        for(int i = 0; i < PACK_SIZE; ++i){
            if(data[i]!= NULL){
                delete data[i];
                data[i] = NULL;
            }
        }
        delete[] data;
        data = NULL;
    }
}

void SingleEndProcessor::consumerTask(ThreadConfig* config){
    while(true) {
        if(config->canBeStopped()){
            ++mFinishedThreads;
            break;
        }
        while(mRepo.writePos <= mRepo.readPos) {
            if(mProduceFinished)
                break;
            usleep(1000);
        }

        if(mProduceFinished && mRepo.writePos == mRepo.readPos){
            ++mFinishedThreads;
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " data processing completed";
                loginfo(msg);
            }
            break;
        }
        
        if(mProduceFinished){
            if(mOptions->verbose) {
                string msg = "thread " + to_string(config->getThreadId() + 1) + " is processing the " + to_string(mRepo.readPos) + " / " + to_string(mRepo.writePos) + " pack";
                loginfo(msg);
            }
            consumePack(config);
        } else {
            consumePack(config);
        }
    }

    if(mFinishedThreads == mOptions->thread) {
        if(mLeftWriter)
            mLeftWriter->setInputCompleted();
        if(mFailedWriter)
            mFailedWriter->setInputCompleted();
    }

    if(mOptions->verbose) {
        string msg = "thread " + to_string(config->getThreadId() + 1) + " finished";
        loginfo(msg);
    }
}

void SingleEndProcessor::writeTask(WriterThread* config)
{
    while(true) {
        if(config->isCompleted()){
            // last check for possible threading related issue
            config->output();
            break;
        }
        config->output();
    }

    if(mOptions->verbose) {
        string msg = config->getFilename() + " writer finished";
        loginfo(msg);
    }
}
