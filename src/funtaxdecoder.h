#ifndef FUNTAXDECODER_H
#define FUNTAXDECODER_H

#include <zlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <vector>

#include "common.h"
#include "phylotree.h"
#include "funtaxoptions.h"

extern std::mutex mtxTreR;
extern std::mutex mtxTreW;

using namespace std;

struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        size_t seed = 0;
        for (const auto& elem : v) {
            seed ^= std::hash<int>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
class FunTaxDecoder{
    public:
        FunTaxDecoder(PhyloOptions *& mOptions);
        ~FunTaxDecoder();

    public:
       void process();

    private:
        void readFunTax();
        void decode();
        void decodeEach();
        std::pair<std::string, std::string> decodeFunTax(std::unordered_set<std::string> &locSet);
        std::string decodeTax(std::unordered_set<std::string>& locSet);
        std::string decodeFun(std::unordered_set<std::string>& locSet);
        void decodeTaxonSample(std::map<std::string, std::map<std::string, uint32>>& tTaxMap);
        void decodeFunSample(std::map<std::string, std::map<std::string, uint32>>& tFunMap, 
            std::map<std::string, std::map<std::string, uint32>>& tPureFunMap, 
            std::map<std::string, std::map<std::string, uint32>>& tGeneFunMap);

    private:
        PhyloOptions *mOptions;
        const int buffer_size = 8192;
        std::unordered_map<std::string, std::unordered_map<std::string, uint32>> totFTMap;
        std::unordered_map<std::string, uint32> ftMap;
        std::unordered_map<std::string, uint32> samFunMap;
        std::unordered_map<std::string, uint32> samTaxMap;
        std::unordered_set<std::string> ftSet;
        std::queue<std::string> ftQueue;
        std::unordered_set<std::string> samFunSet;
        std::unordered_set<std::string> samTaxSet;
        PhyloTree* mPhyloTree;
        FunTaxFreq* mFunTaxFreq;
        std::unordered_map<std::string, std::pair<std::string, std::string>> mFunTaxPair;//line, tax, fun
        std::set<std::string> uniqFuns;
        std::set<std::string> uniqPureFuns;
        std::set<std::string> uniqGeneFuns;
        std::set<std::string> uniqTaxons;
        std::map<char, uint8_t> taxRankMap;
};

#endif /* FUNTAXDECODER_H */