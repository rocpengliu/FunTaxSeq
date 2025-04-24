#include "funtaxoptions.h"
#include <regex>

 PhyloOptions::PhyloOptions(){
     db = "";
     prefix = "";
     outTaxon = "";
     outFun = "";
     outGeneFun = "";
     outPureFun = "";
     gTree = "";
     tTree = "";
     taxonGenomeSize = "";
     orthGeneSize = "";
     geneAno = "";
     orthAno = "";
     geneDNADup = "";
     geneProDup = "";
     sampleDir = "";
     verbose = false;
     debug = false;
     samples.clear();
     thread = 4;
     outTree = "";
     taxLevels = {"kindom", "phylum", "class", "order", "family", "genus", "species"};
 }

PhyloOptions::~PhyloOptions(){
}

void PhyloOptions::parseSample(){
    if(sampleDir.empty()) error_exit("sample directory is not specified");
    if(!is_directory(sampleDir)) error_exit(sampleDir + " sample directory does not exist");
    DIR* dir = opendir(sampleDir.c_str());
    std::string pattern = "_funtax.txt.gz"; // Change the pattern as needed
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string filename = entry->d_name;
        if (filename != "." && filename != ".." && filename.find(pattern) != std::string::npos) {
            samples.push_back(sampleDir + "/" + filename);
        }
    }
    closedir(dir);
    if(samples.empty()) error_exit("no sample file *_funtax.txt.gz found!");
}

GeneNode::GeneNode(){
    id = "";
    geneSize = 0;
    par = "";
    taxon = "";
    taxonLev = 0;
    anno = "";
    koSet.clear();
    goSet.clear();
}

GeneNode::~GeneNode(){
}

std::string GeneNode::print(std::string idpar, std::string type){
    bool goGo = false;
    bool goKo = false;
    std::stringstream ss;
    //ss << (idpar == "id" ? id : par) << "|" << anno << "|" << taxon;
    ss << (idpar == "id" ? id : par) << "|" << anno;
    if(type == "gene"){

    } else {
        if(type == "go"){
            goGo = true;
        } else if(type == "ko"){
            goKo = true;
        } else {
            goGo = true;
            goKo = true;
        }
    }

    if(goGo){
        ss << "|";
        if(goSet.empty()){
            ss << 0;
        } else {
                for (const auto & it : goSet){
                    if(it == 0){
                        ss << 0;
                    } else {
                        ss << "GO:" << std::setw(7) << std::setfill('0') << it << (it == *goSet.rbegin() ? "" : ";");
                    }
                }
        }
    }

    if(goKo){
        ss << "|";
        if (koSet.empty()){
            ss << 0;
        } else{
            for (const auto &it : koSet){
                if (it == 0){
                    ss << 0;
                } else{
                    ss << "K" << std::setw(5) << std::setfill('0') << it << (it == *koSet.rbegin() ? "" : ";");
                }
            }
        }
    }
    return ss.str();
}

std::string GeneNode::print2(std::string idpar){
    return (idpar == "id" ? id : par);
}

std::string GeneNode::print3(bool printGeneSize = false){
    bool goGo = true;
    bool goKo = true;
    std::stringstream ss;
    std::string pure_anno = std::regex_replace(anno, std::regex("\\["), "(");
    pure_anno = std::regex_replace(pure_anno, std::regex("\\]"), ")");
    ss << pure_anno << "|" << taxon << "|" << id;
    if(goGo){
        ss << "|";
        if(goSet.empty()){
            ss << 0;
        } else {
                for (const auto & it : goSet){
                    if(it == 0){
                        ss << 0;
                    } else {
                        //ss << it << (it == *goSet.rbegin() ? "" : ";");
                        ss << "GO:" << std::setw(7) << std::setfill('0') << it << (it == *goSet.rbegin() ? "" : ";");
                    }
                }
        }
    }

    if(goKo){
        ss << "|";
        if (koSet.empty()){
            ss << 0;
        } else{
            for (const auto &it : koSet){
                if (it == 0){
                    ss << 0;
                } else{
                    //ss << it << (it == *koSet.rbegin() ? "" : ";");
                    ss << "K" << std::setw(5) << std::setfill('0') << it << (it == *koSet.rbegin() ? "" : ";");
                }
            }
        }
    }
    if(printGeneSize)
        ss << "|" << geneSize;
    return ss.str();
}


SimGeneNode::SimGeneNode(){
    id = "";
    anno = "";
    koSet.clear();
    goSet.clear();
}

SimGeneNode::SimGeneNode(std::string id){
    this->id = id;
    anno = "";
    koSet.clear();
    goSet.clear();
}

SimGeneNode::~SimGeneNode(){
}

std::string SimGeneNode::print(){
    std::stringstream ss;
    ss << id << "\t" << anno << "\tGO:";
    for (const auto & it : goSet){
        ss << it << ";";
    }
    ss << "\tKO:";
    for (const auto &it : koSet){
        ss << it << ";";
    }
    return ss.str();
}