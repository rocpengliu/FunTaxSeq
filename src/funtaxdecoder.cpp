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
    // uniqPureFuns.clear();
    // uniqGeneFuns.clear();
    uniqTaxons.clear();
    pureFunSizeCountPairMap.clear();
    geneSizeCountPairMap.clear();
    // taxRankMap = {{'k', 0},{'p', 1},{'c', 2},{'o', 3},{'f', 4},{'g', 5},{'s', 6},{'t', 7}};
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
    int tot_samples = mOptions->samples.size() + 1;
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
                                auto tmp_id_set = splitStrInt2<std::set, std::string>(tmp_id->second);
                                id_set.insert(tmp_id_set.begin(), tmp_id_set.end());
                            }
                        } else if(split_vec[1] == "pro"){
                            split_vec2 = splitStr(split_vec.at(2));
                            id_set.insert(split_vec2.begin(), split_vec2.end());
                            for(const auto & itv : split_vec2){
                                auto tmp_id = mPhyloTree->geneProDupMap.find(itv);
                                if(tmp_id == mPhyloTree->geneProDupMap.end()) continue;
                                auto tmp_id_set = splitStrInt2<std::set, std::string>(tmp_id->second);
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
    std::map<std::string, std::map<std::string, uint32>> tGeneFunMap;
    for(const auto & it : totFTMap){
        for(const auto & it2 : it.second) {
            auto it3 = mFunTaxPair.find(it2.first);
            if(it3 == mFunTaxPair.end())
                continue;
            if(!it3->second.first.empty())
                tTaxMap[it.first][it3->second.first] += it2.second;
            if(!it3->second.second.empty()) {
                tFunMap[it.first][it3->second.second] += it2.second;
                uniqFuns.insert(it3->second.second);

                auto pure_orth_count = getSomeParts(it3->second.second, "gene");
                tPureFunMap[it.first][pure_orth_count.first] += it2.second;
                pureFunSizeCountPairMap[pure_orth_count.first].first += pure_orth_count.second;
                pureFunSizeCountPairMap[pure_orth_count.first].second++;

                auto pure_gene_count = getSomeParts(it3->second.second, "gene_ko_go");
                tGeneFunMap[it.first][pure_gene_count.first] += it2.second;
                geneSizeCountPairMap[pure_gene_count.first].first += pure_gene_count.second;
                geneSizeCountPairMap[pure_gene_count.first].second++;
            }
        }
    }
    std::thread dTaxThread = std::thread(&FunTaxDecoder::decodeTaxonSample, this, std::ref(tTaxMap));
    std::thread dFunThread = std::thread(&FunTaxDecoder::decodeFunSample, this, std::ref(tFunMap), std::ref(tPureFunMap), std::ref(tGeneFunMap));
    if(dTaxThread.joinable()) dTaxThread.join();
    if(dFunThread.joinable()) dFunThread.join();
    if(mOptions->verbose)
        loginfo("decode each sample done!");
}

void FunTaxDecoder::decodeTaxonSample(std::map<std::string, std::map<std::string, uint32>>& tTaxMap){
    std::ofstream *of = new std::ofstream();
    of->open(mOptions->outTaxon.c_str(), std::ofstream::out);
    if(!of->is_open()) error_exit("can not open " + mOptions->outTaxon);
    *of << "#taxon" << "\t" << "genome_size" << "\t";
    for(auto prt = tTaxMap.begin(); prt != tTaxMap.end(); ++prt){
        *of << prt->first << (std::next(prt) == tTaxMap.end() ? "\n" : "\t");
    }
    std::string taxon = "";
    std::unordered_map<std::string, int>::iterator taxon_pair;
    for(const auto & it : uniqTaxons){
        *of << it << "\t";
        taxon.clear();
        taxon = getFirstLastElement(it, false, ';');
        taxon_pair = mPhyloTree->genomeSizeMap.find(taxon);
        *of << (taxon_pair == mPhyloTree->genomeSizeMap.end() ? 0 : taxon_pair->second) << "\t";
        for (auto prt = tTaxMap.begin(); prt != tTaxMap.end(); ++prt){
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
        *of << "#taxon:" << mOptions->taxLevels.at(i) << "\t" << "genome_size" << "\t";
        for(auto sam = subTaxMap.begin(); sam != subTaxMap.end(); ++sam){
            *of << sam->first << (std::next(sam) == subTaxMap.end() ? "\n" : "\t");
        }
        for(const auto & it : subUniqTaxon) {
            *of << it << "\t";
            taxon.clear();
            taxon = getFirstLastElement(it, false, ';');
            taxon_pair = mPhyloTree->genomeSizeMap.find(taxon);
            *of << (taxon_pair == mPhyloTree->genomeSizeMap.end() ? 0 : taxon_pair->second) << "\t";
            for(auto sam = subTaxMap.begin(); sam != subTaxMap.end(); ++sam){
                auto sam2 = sam->second.find(it);
                *of << (sam2 == sam->second.end() ? 0 : sam2->second) << (std::next(sam) == subTaxMap.end() ? "\n" : "\t");
            }
        }
        of->clear();
        of->close();
        if(of){
            delete of;
            of = nullptr;
        }
    }
}

void FunTaxDecoder::decodeFunSample(std::map<std::string, std::map<std::string, uint32>>& tFunMap, 
    std::map<std::string, std::map<std::string, uint32>>& tPureFunMap, 
    std::map<std::string, std::map<std::string, uint32>>& tGeneFunMap){
    std::ofstream* otf = new std::ofstream();
    otf->open(mOptions->outFun.c_str(), std::ofstream::out);
    if(!otf->is_open()) error_exit("can not open " + mOptions->outFun);
    std::string gene = "";
    std::unordered_map<std::string, GeneNode *>::iterator gene_itr;
    *otf << "#ortholog" << "\t" << "gene_size" << "\t";
    for(auto prt = tFunMap.begin(); prt != tFunMap.end(); ++prt){
        *otf << prt->first << (std::next(prt) == tFunMap.end() ? "\n" : "\t");
    }
    for(const auto & it : uniqFuns) {
        *otf << it << "\t";
        // auto itt = mPhyloTree->orthAnoMap.find(it);
        // if(itt == mPhyloTree->orthAnoMap.end()) continue;
        // *otf << itt->second->print3() << "\t";
        // gene.clear();
        // gene = getFirstLastElement(it, true, '|');
        // gene_itr = mPhyloTree->geneAnoMap.find(gene);
        // *otf << it << "\t";
        // bool go = false;
        // if(gene_itr != mPhyloTree->geneAnoMap.end()){
        //     if(gene_itr->second->geneSize != 0){
        //         *otf << gene_itr->second->geneSize << "\t";
        //     } else {
        //         go = true;
        //     }
        // } else {
        //     go = true;
        // }

        // if(go){
        //     gene_itr = mPhyloTree->orthAnoMap.find(gene);
        //     if(gene_itr != mPhyloTree->orthAnoMap.end()){
        //         if(gene_itr->second->geneSize != 0){
        //             *otf << gene_itr->second->geneSize << "\t";
        //         } else {
        //             *otf << 0 << "\t";
        //         }
        //     } else {
        //         *otf << 0 << "\t";
        //     }
        // }

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
    *otf << "#ortholog" << "\t" << "gene_size" << "\t";
    for(auto prt = tPureFunMap.begin(); prt != tPureFunMap.end(); ++prt){
        *otf << prt->first << (std::next(prt) == tPureFunMap.end() ? "\n" : "\t");
    }
    for(const auto & it : pureFunSizeCountPairMap) {
        *otf << it.first << "\t" << (it.second.first / it.second.second) << "\t";
        for (auto pr = tPureFunMap.begin(); pr != tPureFunMap.end(); ++pr){
            auto pr2 = pr->second.find(it.first);
            *otf << (pr2 == pr->second.end() ? 0 : pr2->second) << (std::next(pr) == tPureFunMap.end() ? "\n" : "\t");
        }
    }
    otf->flush();
    otf->clear();
    otf->close();

    otf->open(mOptions->outGeneFun.c_str(), std::ofstream::out);
    if(!otf->is_open()) error_exit("can not open " + mOptions->outGeneFun);
    *otf << "#ortholog" << "\t";
    for(auto prt = tGeneFunMap.begin(); prt != tGeneFunMap.end(); ++prt){
        *otf << prt->first << (std::next(prt) == tGeneFunMap.end() ? "\n" : "\t");
    }
    for(const auto & it : geneSizeCountPairMap) {
        *otf << it.first << "\t" << (it.second.first / it.second.second) << "\t";
        for (auto pr = tGeneFunMap.begin(); pr != tGeneFunMap.end(); ++pr){
            auto pr2 = pr->second.find(it.first);
            *otf << (pr2 == pr->second.end() ? 0 : pr2->second) << (std::next(pr) == tGeneFunMap.end() ? "\n" : "\t");
        }
    }
    otf->flush();
    otf->clear();
    otf->close();

    if(otf){
        delete otf;
        otf = nullptr;
    }
}

std::pair<std::string, std::string> FunTaxDecoder::decodeFunTax(std::unordered_set<std::string>& locSet) {
    std::pair<std::string, std::string> ftp;
    ftp.first = decodeTax(locSet);
    ftp.second = mOptions->gTree.empty() ? decodeFun2(locSet) : decodeFun(locSet);
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
            gene = it->second->print3(true);
        }
        return gene;
    }
    std::unordered_set<std::string> tmpSet;
    for (const auto & it : locSet){
        auto it2 = mPhyloTree->geneAnoMap.find(it);
        if(it2 == mPhyloTree->geneAnoMap.end())
            continue;
        if(it2->second->par == "0"){
            gene_anno_map[it2->second->print3(true)]++;
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
            gene = it2->second->print3(true);
            return gene;
        }
    }
    std::unordered_set<std::string> tmpSet2;
    std::set<tree<std::string*>::iterator, tree<std::string*>::iterator_base_less> treItSet;
    tree<std::string*>::leaf_iterator locf;
    for (const auto & it : tmpSet) {
        locf = std::find_if(mPhyloTree->geneTree->begin_leaf(), mPhyloTree->geneTree->end_leaf(),
                [&it](std::string* & itp) {
                    return *itp == it;
                });
        if (mPhyloTree->geneTree->is_valid(locf)) {
            treItSet.insert(locf);
        } else {
            tmpSet2.insert(it);
        }
    }

    if(!tmpSet2.empty()){
        for (const auto & it : tmpSet2) {
            locf = std::find_if(mPhyloTree->geneTree->begin_post(), mPhyloTree->geneTree->end_post(),
                        [&it](std::string* & itp) {
                            return *itp == it;
                        });
            if (mPhyloTree->geneTree->is_valid(locf)) {
                treItSet.insert(locf);
            }
        }
    }

    //treItSet is size is 1, then just search the ortho map
    if (!treItSet.empty()) {
        //ftNode->funLoc = mPhyloTree->geneTree->lowest_common_ancestor(treItSet);
        auto itt = mPhyloTree->geneTree->lowest_common_ancestor(treItSet);
        auto itt2 = mPhyloTree->orthAnoMap.find(*(itt.node->data));
        if(itt2 != mPhyloTree->orthAnoMap.end()){
            gene = itt2->second->print3(true);
        } else {
            gene = gene2;
        }
        treItSet.clear();
    } else {
        gene = gene2;
    }
    return gene;
}

std::string FunTaxDecoder::decodeFun2(std::unordered_set<std::string>& locSet){
    std::string gene = "";
    if(locSet.size() == 1){
        auto it = mPhyloTree->geneAnoMap.find(*(locSet.begin()));
        if(it != mPhyloTree->geneAnoMap.end()){
            gene = it->second->print3(true);
        }
        return gene;
    }
    std::unordered_map<std::string, int> gene_anno_map;
    std::unordered_set<std::string> parSet;
    for (const auto & it : locSet){
        auto it2 = mPhyloTree->geneAnoMap.find(it);
        if(it2 == mPhyloTree->geneAnoMap.end())
            continue;
        if(it2->second->par == "0"){
            gene_anno_map[it2->second->print3(true)]++;
        } else {
            parSet.insert(it2->second->par);
        }
    }

    auto gene2 = getMapMaxKey(gene_anno_map);
    if(parSet.empty()) {
        return gene2;
    } else if(parSet.size() == 1){
        auto it2 = mPhyloTree->orthAnoMap.find(*(parSet.begin()));
        if(it2 != mPhyloTree->orthAnoMap.end()){
            gene = it2->second->print3(true);
            return gene;
        }
    }
    std::vector<uint8_t> levVec;
    levVec.reserve(parSet.size());
    std::unordered_set<std::string> tmpParSet;
    while(parSet.size() > 1){
        tmpParSet.clear();
        levVec.clear();
        for(const auto & it : parSet){
            auto it2 = mPhyloTree->orthAnoMap.find(it);
            if(it2 == mPhyloTree->orthAnoMap.end()){
                continue;
            }
            levVec.emplace_back(it2->second->taxonLev);
        }
        auto min_lev = std::min_element(levVec.begin(), levVec.end());
        if(min_lev == levVec.end() || static_cast<int>(*min_lev) == 7){
            return "";
        }
        for(const auto & it : parSet){
            auto it2 = mPhyloTree->orthAnoMap.find(it);
            if(it2 == mPhyloTree->orthAnoMap.end()){
                continue;
            }
            if(it2->second->taxonLev == *min_lev){
                if(it2->second->par == "0"){
                    tmpParSet.insert(it2->second->anno);
                } else {
                    tmpParSet.insert(it2->second->par);
                }
            } else {
                tmpParSet.insert(it);
            }
        }
        if(tmpParSet.empty()){
            return "";
        } else if(tmpParSet.size() == 1){
            auto it2 = mPhyloTree->orthAnoMap.find(*(tmpParSet.begin()));
            if(it2 == mPhyloTree->orthAnoMap.end()){
                return "";
            } else {
                return it2->second->print3(true);
            }
        } else {
            std::swap(parSet, tmpParSet);
            tmpParSet.clear();
            levVec.clear();
        }
    }
}