#ifndef PHYLOTREE_H
#define PHYLOTREE_H

#include <zlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <unordered_map>
#include <cstdint> 
#include <algorithm>
#include <future>

#include "util.h"
#include "tree.h"
#include "funtaxoptions.h"

extern std::mutex mtxTreR;
extern std::mutex mtxTreW;

class PhyloTree {
public:
    PhyloTree(PhyloOptions *& mOptions);
    PhyloTree(const PhyloTree& orig);
    virtual ~PhyloTree();

public:
    //tree<std::string*>* taxonTree;
    //tree<std::string*>* geneTree;
    std::unordered_map<std::string, GeneNode*> geneAnoMap;
    std::unordered_map<std::string, GeneNode*> orthAnoMap;
    std::unordered_map<std::string, std::string> geneDNADupMap;
    std::unordered_map<std::string, std::string> geneProDupMap;
    std::unordered_map<std::string, std::pair<std::string, uint16_t>> markerTaxonSizeMap;
    std::unordered_map<std::string, int> genomeSizeMap;
    //std::unordered_map<std::string, int> geneSizeMap;
    std::shared_ptr<tree<std::string*>> taxonTree;
    //std::shared_ptr<tree<std::string*>> markerTree;
    std::shared_ptr<tree<std::string*>> geneTree;
    std::shared_ptr<tree<SimGeneNode*>> geneNodeTree;

private:
    const int buffer_size = 16384;
    void init();
    void printParKid(std::string tre);
    std::queue<std::string> readGZ(std::string & fl);
    void readGZ(std::string &fl, char type);
    void readAnno(std::queue<std::string> &annoQueue, char type);
    //void readOrthAnno(std::queue<std::string>& orthAnnoQueue);
    void readGeneDup(std::queue<std::string>& geneDupQueue, char type);
    void populateGeneTre();
    // tree<std::string*>* buildTreeLoopPtr(std::string* str);
    // tree<std::string*>* buildTreePtr(std::queue<std::string> & linQue, int numThreads);
    std::shared_ptr<tree<std::string *>> buildTreeLoopPtr(std::string *str);
    std::shared_ptr<tree<std::string *>> buildTreePtr(std::string& db, std::string type);
    std::shared_ptr<tree<SimGeneNode *>> buildTreeLoopPtrNode(std::string *str);
    std::shared_ptr<tree<SimGeneNode *>> buildTreePtrNode(std::string& db);
    std::shared_ptr<tree<uint32 *>> buildTreeLoopIntPtr(std::string *str);
    std::shared_ptr<tree<uint32 *>> buildTreeIntPtr(std::queue<std::string> &linQue, int numThreads);
    void print_children_par(std::shared_ptr<tree<std::string *>>& tre, std::string fpath);
private:
    PhyloOptions * mOptions;
};

#endif /* PHYLOTREE */