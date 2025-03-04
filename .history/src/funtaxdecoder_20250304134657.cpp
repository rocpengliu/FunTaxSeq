#include <chrono> 
#include "tree.h"
#include "funtaxdecoder.h"
#include "util.h"

FunTaxDecoder::FunTaxDecoder(PhyloOptions *& mOptions){
    this->mOptions = mOptions;
    samFunMap.clear();
    samTaxMap.clear();
    ftMap.clear();
    ftSet.clear();
    samFunSet.clear();
    samTaxSet.clear();
    totFTMap.clear();
    mPhyloTree = new PhyloTree(mOptions);
    mFunTaxFreq = new FunTaxFreq();
    mFunTaxPair.clear();
    uniqFuns.clear();
    uniqPureFuns.clear();
    uniqTaxons.clear();
    //taxRankMap = {{'k', 0},{'p', 1},{'c', 2},{'o', 3},{'f', 4},{'g', 5},{'s', 6},{'t', 7}};
    taxRankMap = {{'k', 0},{'p', 1},{'c', 2},{'o', 3},{'f', 4},{'g', 5},{'s', 6}};
}

FunTaxDecoder::~FunTaxDecoder(){
    if(mPhyloTree){
        delete mPhyloTree;
        mPhyloTree = NULL;
    }
    if(mFunTaxFreq){
        delete mFunTaxFreq;
        mFunTaxFreq = NULL;
    }
}

void FunTaxDecoder::process(){
    readFunTax();
    if(!totFTMap.empty()){
        loginfo("Number of unique funtax ids: " + std::to_string(ftSet.size()));
        decode();
        decodeEach();
    }
}

void FunTaxDecoder::readFunTax(){
    std::queue<std::string> samQueue;
    for(const auto &it : mOptions->samples){
        samQueue.push(it);
    }
    int numThread = std::min(mOptions->thread, static_cast<int>(mOptions->samples.size()));
    std::thread consumerThreads[numThread];
    if(mOptions->verbose) loginfo("start to read");
    int sample_id = 1;
    int tot_samples = mOptions->samples.size();
    for (int i = 0; i < numThread; ++i){
        consumerThreads[i] = std::thread([this, &samQueue, &sample_id, &tot_samples, &i](){
            std::vector<std::string> split_vec;
            split_vec.reserve(3);
            std::string id = "";
            std::vector<std::string> split_vec2;
            split_vec2.reserve(800);
            std::set<std::string> id_set;
            while(true){
                std::unique_lock<std::mutex> lock(mtxTreR);
                if(samQueue.empty()){
                    lock.unlock();
                    break;
                }
                std::string sample = samQueue.front();
                samQueue.pop();
                lock.unlock();
                char buffer[buffer_size];
                std::string line;
                gzFile file = gzopen(sample.c_str(), "rb");
                if (!file) error_exit("Error: Failed to open file " + sample);
                std::unordered_map<std::string, uint32> tmpMap;
                std::unordered_set<std::string> tmpSet;
                while (gzgets(file, buffer, buffer_size) != NULL){
                    line = buffer; // Output or process each line
                    if (line.size() > 0){
                        trimEnds(&line);
                        split_vec.clear();
                        splitStr(line, split_vec);
                        if(split_vec.size() != 3){
                            continue;
                        }
                        split_vec2.clear();
                        id_set.clear();
                        if(split_vec[1] == "host"){
                            continue;
                        } else if(split_vec[1] == "dna"){
                            split_vec2 = splitStr(split_vec.at(2));
                            id_set.insert(split_vec2.begin(), split_vec2.end());
                            for(const auto & itv : split_vec2){
                                auto tmp_id = mPhyloTree->geneDNADupMap.find(itv);
                                if(tmp_id == mPhyloTree->geneDNADupMap.end()) continue;
                                auto tmp_id_set = splitStrInt<std::set, std::string>(tmp_id->second, ";");
                                id_set.insert(tmp_id_set.begin(), tmp_id_set.end());
                            }
                        } else if(split_vec[1] == "pro"){
                            split_vec2 = splitStr(split_vec.at(2));
                            id_set.insert(split_vec2.begin(), split_vec2.end());
                            for(const auto & itv : split_vec2){
                                auto tmp_id = mPhyloTree->geneProDupMap.find(itv);
                                if(tmp_id == mPhyloTree->geneProDupMap.end()) continue;
                                auto tmp_id_set = splitStrInt<std::set, std::string>(tmp_id->second, ";");
                                id_set.insert(tmp_id_set.begin(), tmp_id_set.end());
                            }
                        }
                        id.clear();
                        for(auto its = id_set.begin(); its != id_set.end(); ++its){
                            id += (*its + ";");
                        }
                        tmpMap[id]++;
                        tmpSet.insert(id);
                    }
                }
                gzclose(file);
                std::unique_lock<std::mutex> lock2(mtxTreW);
                totFTMap[removeExtension(basename(sample), "_funtax.txt.gz")] = tmpMap;
                ftSet.insert(tmpSet.begin(), tmpSet.end());
                ++sample_id;
                lock2.unlock();
                loginfo("File " + std::to_string(sample_id) + "/" + std::to_string(tot_samples) + " " + removeExtension(basename(sample), "_funtax.txt.gz") + " readed");
            }
        });
    }
    for(int i = 0; i < numThread; ++i){
        if(consumerThreads[i].joinable()){
            consumerThreads[i].join();
        }
    }
}

void FunTaxDecoder::decode(){
    int numThreads = std::min<int>(mOptions->thread, ftSet.size());
    std::thread consumerThreads[numThreads];
    for (const auto &itm : ftSet){
        ftQueue.push(itm);
    }
    loginfo("Start to decode the funtax ids");
    auto startTime = std::chrono::high_resolution_clock::now();
    uint32 numIds = ftQueue.size();
    for(int i = 0; i < numThreads; ++i){
        consumerThreads[i] = std::thread([this, &i, &numIds, &startTime](){
            std::unordered_set<std::string> locSet;
            std::unordered_map<std::string, std::pair<std::string, std::string>> mFunTaxPairSub;
            while(true){
                std::unique_lock<std::mutex> lock(mtxTreR);
                if (ftQueue.empty()) {
                    lock.unlock();
                    break;
                }
                std::string ft = ftQueue.front();
                ftQueue.pop();
                size_t decoded = numIds - ftQueue.size();
                lock.unlock();
                if (decoded % 1000 == 0){
                    auto now = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> elapsed = now - startTime;
                    double progress = static_cast<double>(decoded) / numIds;
                    double timeElapsed = elapsed.count();
                    double estimatedTotalTime = timeElapsed / progress;
                    double remainingTime = estimatedTotalTime - timeElapsed;
                    int remainingHours = static_cast<int>(remainingTime / 3600);
                    int remainingMinutes = static_cast<int>((remainingTime - remainingHours * 3600) / 60);
                    int remainingSeconds = static_cast<int>(remainingTime) % 60;
                    loginfo("decoded " + std::to_string(decoded / 1000) + "k (" +
                            std::to_string(getPer(decoded, numIds)) + "%) funtax ids. Estimated remaining time: " +
                            std::to_string(remainingHours) + "h " +
                            std::to_string(remainingMinutes) + "m " +
                            std::to_string(remainingSeconds) + "s");
                }
                if(ft.empty())
                    continue;
                locSet.clear();
                locSet = splitStr2(ft);
                if(locSet.empty())
                    continue;
                auto pr = decodeFunTax(locSet);
                if(pr.first.empty() && pr.second.empty())
                    continue;
                mFunTaxPairSub[ft] = pr;
            }
            std::unique_lock<std::mutex> lock2(mtxTreW);
            for(const auto & pair : mFunTaxPairSub){
                mFunTaxPair[pair.first] = pair.second;
                if(!pair.second.first.empty()){
                    uniqTaxons.insert(pair.second.first);
                }
                if(!pair.second.second.empty()){
                    uniqFuns.insert(pair.second.second);
                }
            }
            lock2.unlock();
        });
    }

    for (int i = 0; i < numThreads; ++i) {
        if (consumerThreads[i].joinable()) {
            consumerThreads[i].join();
        }
    }
    loginfo("Finished decoding: " + std::to_string(mFunTaxPair.size()));
}

void FunTaxDecoder::decodeEach(){
    if(mOptions->verbose)
        loginfo("start to decode each sample");
    std::map<std::string, std::map<std::string, uint32>> tTaxMap; // sample, tax, count;
    std::map<std::string, std::map<std::string, uint32>> tFunMap; // sample, fun with taxon and orth, count;
    std::map<std::string, std::map<std::string, uint32>> tPureFunMap; // sample, only fun, count;
    for(const auto & it : totFTMap){
        for(const auto & it2 : it.second) {
            auto it3 = mFunTaxPair.find(it2.first);
            if(it3 == mFunTaxPair.end())
                continue;
            if(!it3->second.first.empty()) 
                tTaxMap[it.first][it3->second.first] += it2.second;
            if(!it3->second.second.empty()) {
                tFunMap[it.first][it3->second.second] += it2.second;
                std::string pure_orth = removeNMpart(it3->second.second, 1, 2, '|');
                tPureFunMap[it.first][pure_orth] += it2.second;
                uniqPureFuns.insert(pure_orth);
            }
        }
    }
    std::thread dTaxThread = std::thread(&FunTaxDecoder::decodeTaxonSample, this, std::ref(tTaxMap));
    std::thread dFunThread = std::thread(&FunTaxDecoder::decodeFunSample, this, std::ref(tFunMap), std::ref(tPureFunMap));
    if(dTaxThread.joinable()) dTaxThread.join();
    if(dFunThread.joinable()) dFunThread.join();
    if(mOptions->verbose)
        loginfo("decode each sample done!");
}

void FunTaxDecoder::decodeTaxonSample(std::map<std::string, std::map<std::string, uint32>>& tTaxMap){
    std::ofstream *of = new std::ofstream();
    of->open(mOptions->outTaxon.c_str(), std::ofstream::out);
    if(!of->is_open()) error_exit("can not open " + mOptions->outTaxon);
    *of << "#taxon" << "\t";
    for(auto prt = tTaxMap.begin(); prt != tTaxMap.end(); ++prt){
        *of << prt->first << (std::next(prt) == tTaxMap.end() ? "\n" : "\t");
    }
    for(const auto & it : uniqTaxons){
        *of << it << "\t";
        for(auto prt = tTaxMap.begin(); prt != tTaxMap.end(); ++prt){
            auto prt2 = prt->second.find(it);
            *of << (prt2 == prt->second.end() ? 0 : prt2->second) << (std::next(prt) == tTaxMap.end() ? "\n" : "\t");
        }
    }
    of->clear();
    of->close();
    if(of){
        delete of;
        of = nullptr;
    }

    for(int i = 0; i < mOptions->taxLevels.size(); ++i){
        std::set<string> subUniqTaxon;
        std::map<std::string, std::map<std::string, uint32>> subTaxMap;
        for (const auto & it : tTaxMap){
            for(const auto & it2 : it.second){
                std::string tax = getFirstNsSeps(it2.first, (i + 1));
                if(tax.empty()) continue;
                subTaxMap[it.first][tax] += it2.second;
                subUniqTaxon.insert(tax);
            }
        }
        std::ofstream *of = new std::ofstream();
        of->open(mOptions->prefix + "_taxon_abundance_" + mOptions->taxLevels.at(i) + ".txt", std::ofstream::out);
        if(!of->is_open()) error_exit("can not open " + mOptions->prefix + "_taxon_abundance_" + mOptions->taxLevels.at(i) + ".txt");
        *of << "#taxon:" << mOptions->taxLevels.at(i) << "\t";
        for(auto sam = subTaxMap.begin(); sam != subTaxMap.end(); ++sam){
            *of << sam->first << (std::next(sam) == subTaxMap.end() ? "\n" : "\t");
        }
        for(const auto & it : subUniqTaxon) {
            *of << it << "\t";
            for(auto sam = subTaxMap.begin(); sam != subTaxMap.end(); ++sam){
                auto sam2 = sam->second.find(it);
                *of << (sam2 == sam->second.end() ? 0 : sam2->second) << (std::next(sam) == subTaxMap.end() ? "\n" : "\t");
            }
        }
        of->close();
        if(of){
            delete of;
            of = nullptr;
        }
    }
}

void FunTaxDecoder::decodeFunSample(std::map<std::string, std::map<std::string, uint32>>& tFunMap, std::map<std::string, std::map<std::string, uint32>>& tPureFunMap){
    std::ofstream* otf = new std::ofstream();
    otf->open(mOptions->outFun.c_str(), std::ofstream::out);
    if(!otf->is_open()) error_exit("can not open " + mOptions->outFun);
    *otf << "#ortholog" << "\t";
    for(auto prt = tFunMap.begin(); prt != tFunMap.end(); ++prt){
        *otf << prt->first << (std::next(prt) == tFunMap.end() ? "\n" : "\t");
    }
    for(const auto & it : uniqFuns) {
        // auto itt = mPhyloTree->orthAnoMap.find(it);
        // if(itt == mPhyloTree->orthAnoMap.end()) continue;
        // *otf << itt->second->print3() << "\t";
        *otf << it << "\t";
        for (auto pr = tFunMap.begin(); pr != tFunMap.end(); ++pr){
            auto pr2 = pr->second.find(it);
            *otf << (pr2 == pr->second.end() ? 0 : pr2->second) << (std::next(pr) == tFunMap.end() ? "\n" : "\t");
        }
    }
    otf->flush();
    otf->clear();
    otf->close();

    otf->open(mOptions->outPureFun.c_str(), std::ofstream::out);
    if(!otf->is_open()) error_exit("can not open " + mOptions->outPureFun);
    *otf << "#ortholog" << "\t";
    for(auto prt = tPureFunMap.begin(); prt != tPureFunMap.end(); ++prt){
        *otf << prt->first << (std::next(prt) == tPureFunMap.end() ? "\n" : "\t");
    }
    for(const auto & it : uniqPureFuns) {
        *otf << it << "\t";
        for (auto pr = tPureFunMap.begin(); pr != tPureFunMap.end(); ++pr){
            auto pr2 = pr->second.find(it);
            *otf << (pr2 == pr->second.end() ? 0 : pr2->second) << (std::next(pr) == tPureFunMap.end() ? "\n" : "\t");
        }
    }
    otf->flush();
    otf->clear();
    otf->close();
    if(otf){
        delete otf;
        otf = nullptr;
    }

    /*

    std::map<std::string, std::map<std::string, uint32>> finalFunMap;
    for(const auto & it : tFunMap){
        std::multimap<std::string, std::pair<std::vector<std::string>, uint32>> annoFunMap;
        std::set<std::string> uniqFunSet;
        for(const auto & it2 : it.second){
            std::string tmpStr = it2.first;
            auto vec = splitStr(tmpStr, "|");
            annoFunMap.insert(std::make_pair(vec.at(0), std::make_pair(vec, it2.second)));
            uniqFunSet.insert(vec.at(0));
        }
        for(const auto & it2 : uniqFunSet){
            auto range = annoFunMap.equal_range(it2);
            uint32 numReadsTmp = 0;
            std::multimap<uint8, std::vector<std::string>> tmpMap;
            for(auto itr = range.first; itr != range.second; ++itr){
                numReadsTmp += itr->second.second;
                if(!itr->second.first.at(1).empty()){
                    auto va = taxRankMap[itr->second.first.at(1)[0]];
                    tmpMap.insert(std::make_pair(va, itr->second.first));
                }
            }
            if(!tmpMap.empty()){
                uint8 rankKey = tmpMap.begin()->first;
                auto rangeIt = tmpMap.equal_range(rankKey);
                int len = std::distance(rangeIt.first, rangeIt.second);
                if(len == 1){
                    finalFunMap[it.first][getStrVec(rangeIt.first->second)] = numReadsTmp;
                } else {
                    std::vector<std::pair<uint32, uint32>> freqVec;
                    uint32 ii = 0;
                    for(auto rIt = rangeIt.first; rIt != rangeIt.second; ++rIt){
                        freqVec.push_back(std::make_pair(ii, static_cast<uint32>(rIt->second.at(2).length() + rIt->second.at(3).length())));
                    }
                    std::sort(freqVec.begin(), freqVec.end(), [](const std::pair<uint32, uint32>& a, const std::pair<uint32, uint32>& b) {
                        return a.second > b.second;
                    });
                    std::advance(rangeIt.first, freqVec.at(0).first);
                    if ( rangeIt.first != rangeIt.second){
                        finalFunMap[it.first][getStrVec(rangeIt.first->second)] = numReadsTmp;
                    }
                }
            }
        }
    */
}

std::pair<std::string, std::string> FunTaxDecoder::decodeFunTax(std::unordered_set<std::string>& locSet) {
    std::pair<std::string, std::string> ftp;
    ftp.first = decodeTax(locSet);
    ftp.second = decodeFun(locSet);
    return ftp;
}

std::string FunTaxDecoder::decodeTax(std::unordered_set<std::string>& locSet) {
    std::set<tree<std::string*>::iterator, tree<std::string*>::iterator_base_less> treItSet;
    tree<std::string*>::leaf_iterator locf;
    std::set<std::string> uniq_taxon_id_set;
    for (const auto & it : locSet) {
        std::string strChar(it);
        size_t pos = strChar.find_first_of(':');
        if (pos == std::string::npos) continue;
        strChar.erase(pos, std::string::npos);
        uniq_taxon_id_set.insert(strChar);
    }
    for(const auto & it : uniq_taxon_id_set){
        locf = std::find_if(mPhyloTree->taxonTree->begin_leaf(),
            mPhyloTree->taxonTree->end_leaf(),
                [&it](std::string* & itp) {
                    return *itp == it;
                });
        if (mPhyloTree->taxonTree->is_valid(locf)) {
            treItSet.insert(locf);
        }
    }
    std::string taxon = "";
    if (!treItSet.empty()) {
        //ftNode->taxLoc = mPhyloTree->taxonTree->lowest_common_ancestor(treItSet);
        taxon = mPhyloTree->taxonTree->lowest_common_ancestor_str(treItSet);
        trimLeft(taxon, "root;");
        treItSet.clear();
    }
    return taxon;
}

std::string FunTaxDecoder::decodeFun(std::unordered_set<std::string>& locSet) {
    std::unordered_map<std::string, int> gene_anno_map;
    std::string gene = "";
    if(locSet.size() == 1){
        auto it = mPhyloTree->geneAnoMap.find(*(locSet.begin()));
        if(it != mPhyloTree->geneAnoMap.end()){
            gene = it->second->print3();
        }
        return gene;
    }
    std::unordered_set<std::string> tmpSet;
    for (const auto & it : locSet){
        auto it2 = mPhyloTree->geneAnoMap.find(it);
        if(it2 == mPhyloTree->geneAnoMap.end())
            continue;
        if(it2->second->par == "0"){
            gene_anno_map[it2->second->print3()]++;
        } else {
            tmpSet.insert(it2->second->par);
        }
    }

    auto gene2 = getMapMaxKey(gene_anno_map);
    if(tmpSet.empty()) {
        return gene2;
    } else if(tmpSet.size() == 1){
        auto it2 = mPhyloTree->orthAnoMap.find(*(tmpSet.begin()));
        if(it2 != mPhyloTree->orthAnoMap.end()){
            gene = it2->second->print3();
            return gene;
        }
    }
    std::set<tree<std::string*>::iterator, tree<std::string*>::iterator_base_less> treItPathSet;
    std::set<tree<std::string*>::iterator, tree<std::string*>::iterator_base_less> treItSet;
    tree<std::string*>::leaf_iterator locf;
    for (const auto & it : tmpSet) {
        locf = std::find_if(mPhyloTree->geneTree->begin_leaf(),
                mPhyloTree->geneTree->end_leaf(),
                [&it](std::string* & itp) {
                    return *itp == it;
                });
        if (mPhyloTree->geneTree->is_valid(locf)) {
            treItSet.insert(locf);
        } else {
            auto path = mPhyloTree->geneTree->path_from_iterator(locf, mPhyloTree->geneTree->begin());
            path.erase(path.begin() + 2, path.end());
            tree<std::string *>::iterator itPath = mPhyloTree->geneTree->iterator_from_path(path, mPhyloTree->geneTree->begin());
            treItPathSet.insert(itPath);
        }
    }

    if(!tmpSet2.empty()){
        if(treItSet.empty()){
            tree<std::string*>::post_order_iterator locfp;
            for (const auto & it : tmpSet2) {
                locf = std::find_if(mPhyloTree->geneTree->begin_post(), mPhyloTree->geneTree->end_post(),
                        [&it](std::string* & itp) {
                            return *itp == it;
                        });
                if (mPhyloTree->geneTree->is_valid(locf)) {
                    treItSet.insert(locf);
                }
            }
        } else {
            for(const auto & itp : treItSet){

            }
        }
    }

    // if(treItSet.empty()){
    //     tree<std::string*>::pre_order_iterator locf2;
    //     for (const auto & it : tmpSet) {
    //         locf2 = std::find_if(mPhyloTree->geneTree->begin(), mPhyloTree->geneTree->end(),
    //                 [&it](std::string* & itp) {
    //                     return *itp == it;
    //                 });
    //         if (mPhyloTree->geneTree->is_valid(locf2)) {
    //             treItSet.insert(locf2);
    //         }
    //     }
    // }

    //treItSet is size is 1, then just search the ortho map
    if (!treItSet.empty()) {
        //ftNode->funLoc = mPhyloTree->geneTree->lowest_common_ancestor(treItSet);
        auto itt = mPhyloTree->geneTree->lowest_common_ancestor(treItSet);
        auto itt2 = mPhyloTree->orthAnoMap.find(*(itt.node->data));
        if(itt2 != mPhyloTree->orthAnoMap.end()){
            //gene = itt2->second->id;
            //gene = itt2->first;
            gene = itt2->second->print3();
            // else {
            //     gene = gene2;
            // }
        } else {
            gene = gene2;
        }
        treItSet.clear();
    } else {
        gene = gene2;
    }
    return gene;
}

/*
std::string FunTaxDecoder::decodeFun(std::unordered_set<std::string>& locSet) {
    std::string gene = "";
    if(locSet.size() == 1){
        auto it = mPhyloTree->geneAnoMap.find(*(locSet.begin()));
        if(it != mPhyloTree->geneAnoMap.end()){
            gene = it->second->print2("par");
        }
        return gene;
    }
    std::unordered_set<std::string> tmpSet;
    for (const auto & it : locSet){
        auto it2 = mPhyloTree->geneAnoMap.find(it);
        if(it2 == mPhyloTree->geneAnoMap.end())
            continue;
        if(it2->second->par != "0") tmpSet.insert(it2->second->par);
    }

    std::set<tree<std::string*>::iterator, tree<std::string*>::iterator_base_less> treItSet;
    tree<std::string*>::leaf_iterator locf;
    for (const auto & it : tmpSet) {
        locf = std::find_if(mPhyloTree->geneTree->begin_leaf(),
                mPhyloTree->geneTree->end_leaf(),
                [&it](std::string* & itp) {
                    return *itp == it;
                });
        if (mPhyloTree->geneTree->is_valid(locf)) {
            treItSet.insert(locf);
        }
    }
    
    if (!treItSet.empty()) {
        //ftNode->funLoc = mPhyloTree->geneTree->lowest_common_ancestor(treItSet);
        auto itt = mPhyloTree->geneTree->lowest_common_ancestor(treItSet);
        auto itt2 = mPhyloTree->orthAnoMap.find(*(itt.node->data));
        if(itt2 != mPhyloTree->orthAnoMap.end()){
            gene = itt2->second->print2();
            //gene = itt2->first;
        }
        treItSet.clear();
    }
    return gene;
}
*/