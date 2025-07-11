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
#include <utility>

#include "common.h"
#include "phylotree.h"
#include "funtaxoptions.h"

extern std::mutex mtxTreR;
extern std::mutex mtxTreW;

using namespace std;

// struct VectorHash {
//     size_t operator()(const std::vector<int>& v) const {
//         size_t seed = 0;
//         for (const auto& elem : v) {
//             seed ^= std::hash<int>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//         }
//         return seed;
//     }
// };
class FunTaxDecoder{
    public:
        FunTaxDecoder(PhyloOptions *& mOptions);
        ~FunTaxDecoder();

    public:
       void process();

    private:
        void readFunTax();
        void decode();
        void decodeMarker();
        void decodeEach();
        void decodeEachMarker();
        std::pair<std::string, std::string> decodeFunTax(std::unordered_set<std::string> &locSet);
        std::string decodeTax(std::unordered_set<std::string>& locSet);
        std::string decodeMarkerTax(std::unordered_set<std::string>& locSet);
        std::string decodeFun(std::unordered_set<std::string>& locSet);
        std::string decodeFun2(std::unordered_set<std::string>& locSet);
        void decodeTaxonSample(std::map<std::string, std::map<std::string, uint32>>& tTaxMap);
        void decodeFunSample(std::map<std::string, std::map<std::string, uint32>>& tFunMap,
            std::map<std::string, std::map<std::string, uint32>>& tPureFunMap,
            std::map<std::string, std::map<std::string, uint32>>& tGeneFunMap);

    private:
        PhyloOptions *mOptions;
        const int buffer_size = 8192;
        std::unordered_map<std::string, std::unordered_map<std::string, uint32>> totFTMap;
        std::unordered_map<std::string, uint32> ftMap;
        std::unordered_map<std::string, std::unordered_map<std::string, uint32>> totMarkerMap;
        std::unordered_map<std::string, uint32> markerMap;
        std::unordered_map<std::string, uint32> samFunMap;
        std::unordered_map<std::string, uint32> samTaxMap;
        std::unordered_set<std::string> ftSet;
        std::unordered_set<std::string> markerSet;
        std::queue<std::string> ftQueue;
        std::queue<std::string> markerQueue;
        std::unordered_set<std::string> samFunSet;
        std::unordered_set<std::string> samTaxSet;
        PhyloTree* mPhyloTree;
        FunTaxFreq* mFunTaxFreq;
        std::unordered_map<std::string, std::pair<std::string, std::string>> mFunTaxPair;//line, tax, fun
        std::unordered_map<std::string, std::string> mMarkerTax;//for marker gene
        std::set<std::string> uniqFuns;
        //std::set<std::string> uniqPureFuns;
        std::map<std::string, std::pair<int, int>> pureFunSizeCountPairMap; //fun with taxon and orth, sum of gene size and number of genes
        //std::set<std::string> uniqGeneFuns;
        std::map<std::string, std::pair<int, int>> geneSizeCountPairMap; //fun with taxon and orth, sum of gene size and number of genes
        std::set<std::string> uniqTaxons;
        std::set<std::string> uniqMarkerTaxons;
};

#endif /* FUNTAXDECODER_H */