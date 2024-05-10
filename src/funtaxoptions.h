#ifndef FUNTAXOPTIONS_H
#define FUNTAXOPTIONS_H

#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <set>
#include <unordered_map>
#include <cstdint> 
#include <dirent.h>

#include "util.h"
#include "tree.h"

using namespace std;
class PhyloOptions{
    public:
        PhyloOptions();
        virtual ~PhyloOptions();

    public:
        void parseSample();

    public:
        std::string db;
        std::string prefix;
        std::string outTaxon;
        std::string outFun;
        std::string gTree;
        std::string tTree;
        std::string geneAno;
        std::string orthAno;
        std::string sampleDir;
        std::string outTree;
        bool verbose;
        bool debug;
        std::vector<std::string> samples;
        int thread;
        std::vector<std::string> taxLevels;
};

class FunTaxFreq{
    public:
        FunTaxFreq(){
            geneFreq.clear();
            taxonFreq.clear();
        }
        void clear(){
            geneFreq.clear();
            taxonFreq.clear();
        }

    public:
        std::unordered_map<std::string*, int> geneFreq;
        std::unordered_map<std::string*, int> taxonFreq;
};

class FunTaxNode{
    public:
    FunTaxNode(){
    }

    public:
        tree<std::string*>::leaf_iterator taxLoc;
        tree<std::string *>::leaf_iterator funLoc;
};

class GeneNode{
    public:
        GeneNode();
        ~GeneNode();

        std::string print(std::string type = "full");

    public:
        std::string id;
        std::string par;
        std::string taxon;
        std::string anno;
        std::set<uint16> koSet;
        std::set<uint32> goSet;
};

class SimGeneNode{
    public:
        SimGeneNode();
        SimGeneNode(std::string id);
        ~SimGeneNode();

        std::string print();

    public:
        std::string id;
        std::string anno;
        std::set<uint16> koSet;
        std::set<uint32> goSet;
};

#endif /* FUNTAXOPTIONS_H */