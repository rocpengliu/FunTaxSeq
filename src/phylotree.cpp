#include "phylotree.h"

PhyloTree::PhyloTree(PhyloOptions *& mOptions) {
    this->mOptions = mOptions;
    geneTree = nullptr;
    geneNodeTree = nullptr;
    taxonTree = nullptr;
    //markerTree = nullptr;
    geneAnoMap.clear();
    orthAnoMap.clear();
    geneDNADupMap.clear();
    geneProDupMap.clear();
    genomeSizeMap.clear();
    //markerSizeMap.clear();
    markerTaxonSizeMap.clear();
    // geneSizeMap.clear();
    init();
}

PhyloTree::PhyloTree(const PhyloTree& orig) {
}

PhyloTree::~PhyloTree() {
    std::thread dTTreThread = std::thread([this](){
        if (mOptions->verbose)
            cerr << "start to free taxon tree" << "\n";
        if(taxonTree){
            for (auto it = taxonTree->begin(); it != taxonTree->end(); ++it){
                if (it.node->data){
                    delete it.node->data;
                    it.node->data = nullptr;
                }
            }
        }
        if(mOptions->verbose) cerr << "taxon tree free finished" << "\n";
    });

    // std::thread dMTreThread = std::thread([this](){
    //     if (mOptions->verbose)
    //         cerr << "start to free marker tree" << "\n";
    //     if(markerTree){
    //         for (auto it = markerTree->begin(); it != markerTree->end(); ++it){
    //             if (it.node->data){
    //                 delete it.node->data;
    //                 it.node->data = nullptr;
    //             }
    //         }
    //     }
    //     if(mOptions->verbose) cerr << "marker tree free finished" << "\n";
    // });

    std::thread dFTreThread;
    if(!mOptions->gTree.empty()){
        dFTreThread = std::thread([this](){
            if(mOptions->verbose) cerr << "start to free gene ortholog tree" << "\n";
            for(auto it = geneTree->begin(); it != geneTree->end(); ++it){
                if(it.node->data){
                delete it.node->data;
                it.node->data = nullptr;
                }
            }
            //delete geneTree;
            if(mOptions->verbose) cerr << "gene ortholog tree free finished" << "\n";
        });
    }

/*
    std::thread dGenThread = std::thread([this](){
        if(mOptions->verbose) cerr << "start to free gene annotation map" << "\n";
        int ii = 0;
        //do not perform this, as this takes too long time, and won't cause problem as it will be deleted
        for(auto it = geneAnoMap.begin(); it != geneAnoMap.end(); ++it){
            ++ii;
            if(it->second){
                if(ii % 1000000 == 0) loginfo("deleted gene map: ", ii);
                delete it->second;
            //it->second = nullptr;
            }
        }
        geneAnoMap.clear();
        if(mOptions->verbose) cerr << "gene annotation map free finished!" << "\n";
    });
    */

    std::thread dOrthThread = std::thread([this](){
        if(mOptions->verbose) cerr << "start to free ortho annotation map" << "\n";//could remove as this is already deleted.
        for(auto it = orthAnoMap.begin(); it != orthAnoMap.end(); ++it){
            if(it->second){
            delete it->second;
            //it->second = nullptr;
            }
        }
        orthAnoMap.clear();
        if(mOptions->verbose) cerr << "ortho annotation map free finished!" << "\n";
    });

    const int chunkSize = geneAnoMap.size() / mOptions->thread;
    std::vector<std::future<void>> futures;
    futures.reserve(mOptions->thread);
    auto it = geneAnoMap.begin();

    for (int i = 0; i < mOptions->thread; ++i) {
        auto startIt = it;
        std::advance(it, chunkSize);
        if (i == mOptions->thread - 1) {
            it = geneAnoMap.end();
        }
        futures.push_back(std::async(std::launch::async, [startIt, it] {
            for (auto iter = startIt; iter != it; ++iter) {
                if (iter->second != nullptr) {
                    delete iter->second;
                }
            }
        }));
    }

    for (auto& future : futures) {
        future.get();
    }
    geneAnoMap.clear();

    if(dTTreThread.joinable()){
        dTTreThread.join();
    }
    // if(dMTreThread.joinable()){
    //     dMTreThread.join();
    // }
    if(dFTreThread.joinable()){
        dFTreThread.join();
    }
    //if(dGenThread.joinable()){
        //dGenThread.join();
    //}
    if(dOrthThread.joinable()){
        dOrthThread.join();
    }
    loginfo("finished all the deleting");
}

void PhyloTree::init(){
    std::queue<std::string> geneAnnoQueue;
    std::thread readGenoThread = std::thread([this, &geneAnnoQueue]{
        geneAnnoQueue = readGZ(mOptions->geneAno);
    });

    std::queue<std::string> orthAnnoQueue;
    std::thread readOrthThread = std::thread([this, &orthAnnoQueue]{
        orthAnnoQueue = readGZ(mOptions->orthAno);
    });

    //std::queue<std::string> geneDNADupQueue;
    std::thread geneDNADupThread = std::thread([this]{
        readGZ(mOptions->geneDNADup, 'd');
    });

    //std::queue<std::string> geneProDupQueue;
    std::thread geneProDupThread = std::thread([this]{
        readGZ(mOptions->geneProDup, 'p');
    });

    std::thread taxonGenomeSizeThread = std::thread([this]{
        readGZ(mOptions->taxonGenomeSize, 't');
    });

    std::thread markerTaxonThread = std::thread([this]{
        readGZ(mOptions->markerTaxonSize, 'm');
    });

    if(readGenoThread.joinable()){
        readGenoThread.join();
    }

    if(readOrthThread.joinable()){
        readOrthThread.join();
    }
    if(geneDNADupThread.joinable()){
        geneDNADupThread.join();
    }
    if(geneProDupThread.joinable()){
        geneProDupThread.join();
    }
    if(taxonGenomeSizeThread.joinable()){
        taxonGenomeSizeThread.join();
    }
    if(markerTaxonThread.joinable()){
        markerTaxonThread.join();
    }
    // if(orthGeneSizeThread.joinable()){
    //     orthGeneSizeThread.join();
    // }
    // readGeneAnno(geneAnnoQueue);
    // readOrthAnno(orthAnnoQueue);
    // readGeneDup(geneDNADupQueue, 'd');
    // readGeneDup(geneProDupQueue, 'p');
    readAnno(geneAnnoQueue, 'g');
    readAnno(orthAnnoQueue, 'o');
    if (mOptions->verbose)
        loginfo("start to build taxon tree!");
    taxonTree = buildTreePtr(mOptions->tTree, "taxon");
    if(mOptions->verbose) loginfo("finished to build taxon tree!");
    if(taxonTree->size() < 3) error_exit("built taxon tree size must be no less than 2: ");
    if(mOptions->verbose) cerr << "taxon tree size is " << taxonTree->size() << " and has " << taxonTree->begin().number_of_descent() << " descents" << "\n";
    //print_children_par(taxonTree, "taxon_tree_par_children.txt");

    if(!mOptions->gTree.empty()){
        if(mOptions->verbose) loginfo("start to build gene ortholog tree!");
        geneTree = buildTreePtr(mOptions->gTree, "ortholog");
        if(mOptions->verbose) loginfo("finished to build gene ortholog tree!");
        if(geneTree->size() < 1) error_exit("built gene tree size must be no less than 1: ");
        if(mOptions->verbose) cerr << "gene ortholog tree size is " << geneTree->size() << " and has " << geneTree->begin().number_of_descent() << " descents" << "\n";
    }

    if (mOptions->verbose)
        loginfo("all the initiation completed!");
    // if(mOptions->marker){
    //     if(mOptions->verbose) loginfo("start to build marker tree!");
    //     markerTree = buildTreePtr(mOptions->mTree, "marker");
    //     if(mOptions->verbose) loginfo("finished to build marker tree!");
    //     if(geneTree->size() < 3) error_exit("built marker tree size must be no less than 3: ");
    //     if(mOptions->verbose) cerr << "gene marker tree size is " << markerTree->size() << " and has " << markerTree->begin().number_of_descent() << " descents" << "\n";
    // }
}

std::shared_ptr<tree<std::string*>> PhyloTree::buildTreeLoopPtr(std::string* str) {
    std::shared_ptr<tree<std::string*>> tre(new tree<std::string*>());
    size_t pos = str->find_first_of("(");
    if (pos == std::string::npos) {
        error_exit("no first '(' have been found!");
    }
    std::string root = str->substr(0, pos);
    str->erase(0, pos);
    tree<std::string*>::leaf_iterator loc = tre->set_head(new std::string(root));
    while (!str->empty()) {
        char leftC = str->at(0); //left char must be one of the (, );
        size_t pos = str->find_first_of("()", 1);
        if (pos == std::string::npos) {
            loc = tre->parent(loc);
            break;
        }
        std::string locStr = str->substr(1, pos - 1);
        char rightC = str->at(pos); //right char must be one of the (, );
        str->erase(0, pos);
        if (leftC == '(') {
            auto locVec = splitStr(locStr);
            if (locVec.empty()) {
                error_exit("must have sth after1 '('!");
            }
            loc = tre->append_child(loc, new std::string(locVec.at(0)));
            for (int s = 1; s != locVec.size(); ++s) {
                loc = tre->insert_after(loc, new std::string(locVec.at(s)));
            }
        } else {
            loc = tre->parent(loc);
            auto locVec = splitStr(locStr);
            if (rightC == '(') {
                if (locVec.empty()) {
                    error_exit("must have loc between ')' and '('!");
                }
                for (auto & lo : locVec) {
                    loc = tre->insert_after(loc, new std::string(lo));
                }
            } else {
                for (auto & lo : locVec) {
                    loc = tre->insert_after(loc, new std::string(lo));
                }
            }
        }
    }
    return tre;
}

std::shared_ptr<tree<std::string*>> PhyloTree::buildTreePtr(std::string & db, std::string type) {
    std::ifstream ifs(db.c_str());
    if(!ifs.is_open()) error_exit("can not open tree file: " + db);
    std::string line = "";
    std::queue<std::string> linQue;
    while(std::getline(ifs, line)){
        trimEnds(&line);
        if(line.empty()) continue;
        linQue.push(line);
    }
    ifs.close();
    if(mOptions->verbose) cerr << linQue.size() << " lines are readed for " << type << " tree" << "\n";

    std::shared_ptr<tree<std::string*>> tre(new tree<std::string*>());
    tre->set_head(new std::string("root"));
    auto locDummy = tre->append_child(tre->begin(), new std::string("dummy"));
    cCout("initiate a super root tree", tre->size());
    int numThreads = std::min(mOptions->thread, static_cast<int>(linQue.size()));
    std::thread consumerThreads[numThreads];
    for (int i = 0; i < numThreads; ++i) {
        consumerThreads[i] = std::thread([this, &i, &linQue, &tre, &locDummy]() {
            while (true) {
                std::unique_lock<std::mutex> lock(mtxTreR);
                if (linQue.empty()) {
                    lock.unlock();
                    break;
                }
                std::string linstr = linQue.front();
                linQue.pop();
                lock.unlock();
                if (linstr.empty() || linstr == "") continue;
                std::shared_ptr<tree<std::string*>> treTmp = buildTreeLoopPtr(&linstr);
                std::unique_lock<std::mutex> lock2(mtxTreW);
                locDummy = tre->insert_subtree(locDummy, treTmp->begin());
                lock2.unlock();
            }
        });
    }
    for (int i = 0; i < numThreads; ++i) {
        if (consumerThreads[i].joinable()) {
            consumerThreads[i].join();
        }
    }
    cCout("final tree depth", tre->max_depth(), 'g');
    return tre;
}

std::shared_ptr<tree<uint32*>> PhyloTree::buildTreeLoopIntPtr(std::string* str) {
    std::shared_ptr<tree<uint32*>> tre(new tree<uint32*>());
    size_t pos = str->find_first_of("(");
    if (pos == std::string::npos) {
        error_exit("no first '(' have been found!");
    }
    uint32 root = static_cast<uint32>(std::stoul(str->substr(0, pos)));
    str->erase(0, pos);
    tree<uint32*>::leaf_iterator loc = tre->set_head(new uint32(root));
    while (!str->empty()) {
        char leftC = str->at(0); //left char must be one of the (, );
        size_t pos = str->find_first_of("()", 1);
        if (pos == std::string::npos) {
            loc = tre->parent(loc);
            break;
        }
        std::string locStr = str->substr(1, pos - 1);
        char rightC = str->at(pos); //right char must be one of the (, );
        str->erase(0, pos);
        if (leftC == '(') {
            auto locVec = splitStr(locStr);
            if (locVec.empty()) {
                error_exit("must have sth after1 '('!");
            }
            loc = tre->append_child(loc, new uint32(static_cast<uint32>(std::stoul((locVec.at(0))))));
            for (int s = 1; s != locVec.size(); ++s) {
                loc = tre->insert_after(loc, new uint32(static_cast<uint32>(std::stoul((locVec.at(s))))));
            }
        } else {
            loc = tre->parent(loc);
            auto locVec = splitStr(locStr);
            if (rightC == '(') {
                if (locVec.empty()) {
                    error_exit("must have loc between ')' and '('!");
                }
                for (auto & lo : locVec) {
                    loc = tre->insert_after(loc, new uint32(static_cast<uint32>(std::stoul(lo))));
                }
            } else {
                for (auto & lo : locVec) {
                    loc = tre->insert_after(loc, new uint32(static_cast<uint32>(std::stoul(lo))));
                }
            }
        }
    }
    return tre;
}

std::shared_ptr<tree<uint32*>> PhyloTree::buildTreeIntPtr(std::queue<std::string> & linQue, int numThreads = 16) {
    cCout("linQue size in the building tree", linQue.size());
    std::shared_ptr<tree<uint32*>> tre(new tree<uint32*>());
    tre->set_head(new uint32(0));
    auto locDummy = tre->append_child(tre->begin(), new uint32(4294967295));
    cCout("initiate a super root tree", tre->size());
    numThreads = std::min(numThreads, static_cast<int> (linQue.size()));
    std::thread consumerThreads[numThreads];
    for (int i = 0; i < numThreads; ++i) {
        consumerThreads[i] = std::thread([this, &i, &linQue, &tre, &locDummy]() {
            while (true) {
                std::unique_lock<std::mutex> lock(mtxTreR);
                if (linQue.empty()) {
                    lock.unlock();
                    break;
                }
                std::string linstr = linQue.front();
                linQue.pop();
                lock.unlock();
                if (linstr.empty() || linstr == "") continue;
                std::shared_ptr<tree<uint32*>> treTmp = buildTreeLoopIntPtr(&linstr);
                std::unique_lock<std::mutex> lock2(mtxTreW);
                locDummy = tre->insert_subtree(locDummy, treTmp->begin());
                lock2.unlock();
            }
        });
    }
    for (int i = 0; i < numThreads; ++i) {
        if (consumerThreads[i].joinable()) {
            consumerThreads[i].join();
        }
    }
    cCout("final tree depth", tre->max_depth(), 'g');
    return tre;
}

/*
tree<std::string*>* PhyloTree::buildTreeLoopPtr(std::string* str) {
    tree<std::string*>* tre = new tree<std::string*>();
    size_t pos = str->find_first_of("(");
    if (pos == std::string::npos) {
        error_exit("no first '(' have been found!");
    }
    std::string root = str->substr(0, pos);
    str->erase(0, pos);
    tree<std::string*>::leaf_iterator loc = tre->set_head(new std::string(root));
    while (!str->empty()) {
        char leftC = str->at(0); //left char must be one of the (, );
        size_t pos = str->find_first_of("()", 1);
        if (pos == std::string::npos) {
            loc = tre->parent(loc);
            break;
        }
        std::string locStr = str->substr(1, pos - 1);
        char rightC = str->at(pos); //right char must be one of the (, );
        str->erase(0, pos);
        if (leftC == '(') {
            auto locVec = splitStr(locStr);
            if (locVec.empty()) {
                error_exit("must have sth after1 '('!");
            }
            loc = tre->append_child(loc, new std::string(locVec.at(0)));
            for (int s = 1; s != locVec.size(); ++s) {
                loc = tre->insert_after(loc, new std::string(locVec.at(s)));
            }
        } else {
            loc = tre->parent(loc);
            auto locVec = splitStr(locStr);
            if (rightC == '(') {
                if (locVec.empty()) {
                    error_exit("must have loc between ')' and '('!");
                }
                for (auto & lo : locVec) {
                    loc = tre->insert_after(loc, new std::string(lo));
                }
            } else {
                for (auto & lo : locVec) {
                    loc = tre->insert_after(loc, new std::string(lo));
                }
            }
        }
    }
    return tre;
}

tree<std::string*>* PhyloTree::buildTreePtr(std::queue<std::string> & linQue, int numThreads = 16) {
    cCout("linQue size in the building tree", linQue.size());
    tree<std::string*>* tre = new tree<std::string*>();
    tre->set_head(new std::string("root"));
    auto locDummy = tre->append_child(tre->begin(), new std::string("dummy"));
    cCout("intitate a super root tree", tre->size());
    numThreads = std::min(numThreads, static_cast<int> (linQue.size()));
    std::thread consumerThreads[numThreads];
    for (int i = 0; i < numThreads; ++i) {
        consumerThreads[i] = std::thread([this, &i, &linQue, &tre, &locDummy]() {
            while (true) {
                std::unique_lock<std::mutex> lock(mtxTreR);
                if (linQue.empty()) {
                    lock.unlock();
                    break;
                }
                std::string linstr = linQue.front();
                linQue.pop();
                lock.unlock();
                if (linstr.empty() || linstr == "") continue;
                tree<std::string*>* treTmp = buildTreeLoopPtr(&linstr);
                std::unique_lock<std::mutex> lock2(mtxTreW);
                locDummy = tre->insert_subtree(locDummy, treTmp->begin());
                lock2.unlock();
            }
        });
    }
    for (int i = 0; i < numThreads; ++i) {
        if (consumerThreads[i].joinable()) {
            consumerThreads[i].join();
        }
    }
    cCout("final tree depth", tre->max_depth(), 'g');
    cCout("final tree size", tre->size(), 'g');
    cCout("final tree children", tre->number_of_children(tre->begin()), 'g');
    cCout("final tree descents", tre->number_of_descent(tre->begin()), 'g');
    return tre;
}
*/
// void PhyloTree::populateGeneTre(){
//     if(mOptions->verbose) loginfo("populate gene tree");
//     tree<SimGeneNode *>::iterator it;
//     for (it = geneNodeTree->begin(); it != geneNodeTree->end(); ++it){
//         auto locn = it.node->data->id;
//         auto loc = orthAnoMap.find(locn);
//         if(loc == orthAnoMap.end())
//             continue;
//         it.node->data->anno = loc->second->anno;
//         it.node->data->goSet = loc->second->goSet;
//         it.node->data->koSet = loc->second->koSet;;
//     }
//     if(mOptions->verbose) loginfo("populate gene tree done");

//     //delete the orthmap to save memeory
// }

void PhyloTree::printParKid(std::string tre){
    if(tre == "taxon"){

    } else {
        std::ofstream *fout = new std::ofstream();
        fout->open(mOptions->outTree.c_str(), std::ofstream::out);
        if(!fout->is_open()){
            error_exit("cannot open file: " + mOptions->outTree);
            delete fout;
            fout = nullptr;
        }
        uint32 n = 0;
        for (auto it = geneTree->begin_breadth_first(); it != geneTree->end_breadth_first(); ++it){
            if(*(it.node->data) == "root") continue;
            auto itp = geneTree->parent(it);
            if(geneTree->is_valid(itp)){
                *fout << *(itp.node->data) << "\t" << *(it.node->data) << "\n";
            }
            ++n;
            if(n % 10000 == 0)
                std::cout << n << std::endl;
        }
        fout->flush();
        fout->clear();
        fout->close();
        if(fout){
            delete fout;
            fout = nullptr;
        }
        cCout("print par kid result finished");
    }
}

void PhyloTree::readAnno(std::queue<std::string>& annoQueue, char type){
    if(mOptions->verbose){
        loginfo("start to populate anno map");
    }
    int numThreads = std::min(static_cast<int>(annoQueue.size()), mOptions->thread);
    std::thread consumerThreads[numThreads];
    for(int i = 0; i < numThreads; ++i){
        consumerThreads[i] = std::thread([this, &annoQueue, &type, &i](){
            while(true){
                std::unique_lock<std::mutex> lock(mtxTreR);
                if(annoQueue.empty()){
                    lock.unlock();
                    break;
                }
                std::string line = annoQueue.front();
                annoQueue.pop();
                if(annoQueue.size() >= 100000 && annoQueue.size() % 100000 == 0){
                    std::string msg = "anno remains: " + std::to_string(annoQueue.size());
                    loginfo(msg, false);
                }
                lock.unlock();
                if(line.empty())
                    continue;
                std::vector<std::string> strVec;
                splitStr(line, strVec);
                if(strVec.size() == 7 && type == 'g'){
                    GeneNode *tmp = new GeneNode();
                    tmp->id = strVec[0];
                    tmp->geneSize = std::stoi(strVec[1]);
                    tmp->par = strVec[2];
                    tmp->taxon = strVec[3];
                    tmp->anno = strVec[4];
                    if(strVec[5] != "0"){
                        tmp->goSet = splitStrInt<std::set, uint32_t>(strVec[5], 'g');
                    }
                    if (strVec[6] != "0"){
                        tmp->koSet = splitStrInt<std::set, uint16_t>(strVec[6], 'k');
                    }
                    std::unique_lock<std::mutex> lock2(mtxTreW);
                    geneAnoMap[strVec[0]] = tmp;
                    lock2.unlock();
                } else if(strVec.size() == 8 && type == 'o'){
                    GeneNode *tmp = new GeneNode();
                    tmp->id = strVec[0];
                    tmp->geneSize = std::stoi(strVec[1]);
                    tmp->par = strVec[2];
                    tmp->taxonLev = static_cast<uint8_t>(strVec[3][0] - '0');
                    tmp->taxon = strVec[4];
                    tmp->anno = strVec[5];
                    if(strVec[6] != "0"){
                        tmp->goSet = splitStrInt<std::set, uint32_t>(strVec[6], 'g');
                    }
                    if (strVec[7] != "0"){
                        tmp->koSet = splitStrInt<std::set, uint16_t>(strVec[7], 'k');
                    }
                    std::unique_lock<std::mutex> lock2(mtxTreW);
                    orthAnoMap[strVec[0]] = tmp;
                    lock2.unlock();
                } else {
                    error_exit(type + " file contain non valid line!");
                }
            }
        });
    }
    for(int i = 0; i < numThreads; i++){
        if(consumerThreads[i].joinable()){
            consumerThreads[i].join();
        }
    }
    if(mOptions->verbose){
        loginfo("read anno map done with size " + std::to_string(type == 'g' ? geneAnoMap.size() : orthAnoMap.size()));
    }
}

// void PhyloTree::readOrthAnno(std::queue<std::string>& orthAnnoQueue){
//     if(mOptions->verbose){
//         loginfo("start to populate ortho map");
//     }
//     int numThreads = std::min(static_cast<int>(orthAnnoQueue.size()), mOptions->thread);
//     std::thread consumerThreads[numThreads];
//     for(int i = 0; i < numThreads; ++i){
//         consumerThreads[i] = std::thread([this, &orthAnnoQueue, &i](){
//             while(true){
//                 std::unique_lock<std::mutex> lock(mtxTreR);
//                 if(orthAnnoQueue.empty()){
//                     lock.unlock();
//                     break;
//                 }
//                 std::string line = orthAnnoQueue.front();
//                 orthAnnoQueue.pop();
//                 if(orthAnnoQueue.size() >= 100000 && orthAnnoQueue.size() % 100000 == 0){
//                     std::string msg = "ortholog remains: " + std::to_string(orthAnnoQueue.size());
//                     loginfo(msg, false);
//                 }
//                 lock.unlock();
//                 if(line.empty())
//                     continue;
//                 std::vector<std::string> strVec;
//                 splitStr(line, strVec);
//                 if(strVec.size() == 6){
//                     GeneNode *tmp = new GeneNode();
//                     tmp->id = strVec[0];
//                     tmp->par = strVec[1];
//                     tmp->taxon = strVec[2];
//                     tmp->anno = strVec[3];
//                     if(strVec[4] != "0"){
//                         tmp->goSet = splitStrInt<std::set, std::string>(strVec[4], ";");
//                     }
//                     if (strVec[5] != "0"){
//                         tmp->koSet = splitStrInt<std::set, std::string>(strVec[5], ";");
//                     }
//                     std::unique_lock<std::mutex> lock2(mtxTreW);
//                     orthAnoMap[strVec[0]] = tmp;
//                     lock2.unlock();
//                 }
//             } });
//     }
//     for(int i = 0; i < numThreads; i++){
//         if(consumerThreads[i].joinable()){
//             consumerThreads[i].join();
//         }
//     }
//     if(mOptions->verbose){
//         loginfo("populate ortho map done with size " + std::to_string(orthAnoMap.size()));
//     }
// }

// void PhyloTree::readGeneDup(std::queue<std::string>& geneDupQueue, char type){
//     if(mOptions->verbose){
//         loginfo("start to populate gene dup map");
//     }
//     int numThreads = std::min(static_cast<int>(geneDupQueue.size()), mOptions->thread);
//     std::thread consumerThreads[numThreads];
//     for(int i = 0; i < numThreads; ++i){
//         consumerThreads[i] = std::thread([this, &geneDupQueue, &i, &type](){
//             while(true){
//                 std::unique_lock<std::mutex> lock(mtxTreR);
//                 if(geneDupQueue.empty()){
//                     lock.unlock();
//                     break;
//                 }
//                 std::string line = geneDupQueue.front();
//                 geneDupQueue.pop();
//                 if(geneDupQueue.size() >= 100000 && geneDupQueue.size() % 100000 == 0){
//                     std::string msg = "gene dup remains: " + std::to_string(geneDupQueue.size());
//                     loginfo(msg, false);
//                 }
//                 lock.unlock();
//                 if(line.empty())
//                     continue;
//                 std::vector<std::string> strVec;
//                 splitStr(line, strVec);
//                 if(strVec.size() == 2){
//                     std::unique_lock<std::mutex> lock2(mtxTreW);
//                     if(type == 'd'){
//                         geneDNADupMap[strVec.at(0)] = strVec.at(1);
//                     } else if(type == 'p') {
//                         geneProDupMap[strVec.at(0)] = strVec.at(1);
//                     }
//                     lock2.unlock();
//                 }
//             } });
//     }
//     for(int i = 0; i < numThreads; i++){
//         if(consumerThreads[i].joinable()){
//             consumerThreads[i].join();
//         }
//     }
//     if(mOptions->verbose){
//         loginfo("populate dup map done with size " + std::to_string(geneDNADupMap.size()));
//     }
// }

std::queue<std::string> PhyloTree::readGZ(std::string & fl){
    std::queue<std::string> lineQueue;
    std::string line = "";
    std::vector<std::string> strVec;
    char buffer[buffer_size];
    gzFile file = gzopen(fl.c_str(), "rb");
    if (!file)
        error_exit("can not open annotation file: " + fl);
    if(mOptions->verbose) loginfo("start to read file " + basename(fl));
    while (gzgets(file, buffer, buffer_size) != NULL){
        line = buffer;
        trimEnds(&line);
        lineQueue.push(line);
    }
    gzclose(file);
    if(mOptions->verbose)
        loginfo("reading file done " + basename(fl) + " " + std::to_string(lineQueue.size()));
    if(lineQueue.empty())
        error_exit(fl + " is empty!");
    return lineQueue;
}

void PhyloTree::readGZ(std::string & fl, char type){
    std::string line = "";
    std::vector<std::string> strVec;
    char buffer[buffer_size];
    gzFile file = gzopen(fl.c_str(), "rb");
    if (!file){
        if(type == 'm' || type == 's'){
            return;
        } else {
            error_exit("can not open annotation file: " + fl);
        }
    }
    if(mOptions->verbose) loginfo("start to read file " + basename(fl));
    uint32_t count = 0;
    while (gzgets(file, buffer, buffer_size) != NULL){
        line = buffer;
        trimEnds(&line);
        strVec.clear();
        splitStr(line, strVec);
        if(strVec.size() == 2){
            if(type == 't'){
                genomeSizeMap[strVec.at(0)] = std::stoi(strVec.at(1));
            } else if(type == 'd'){
                geneDNADupMap[strVec.at(0)] = strVec.at(1);
            } else if(type == 'p'){
                geneProDupMap[strVec.at(0)] = strVec.at(1);
            }
        } else if(strVec.size() == 3){
            if(type == 'm'){
                markerTaxonSizeMap[strVec.at(0)] = std::make_pair(strVec.at(1), static_cast<uint16_t>(std::stoi(strVec.at(2))));
            }
        } else {
            error_exit(basename(fl) + " contains invalid lines!");
        }
        if(count >= 100000 && count % 100000 == 0){
            std::string msg = "read " + formatNumber(count) + " lines for " + type + "!";
            loginfo(msg, false);
        }
        ++count;
    }
    gzclose(file);
    if(mOptions->verbose) loginfo("complete to read file " + basename(fl));
}

std::shared_ptr<tree<SimGeneNode*>> PhyloTree::buildTreeLoopPtrNode(std::string* str) {
    std::shared_ptr<tree<SimGeneNode*>> tre(new tree<SimGeneNode*>());
    size_t pos = str->find_first_of("(");
    if (pos == std::string::npos) {
        error_exit("no first '(' have been found!");
    }
    std::string root = str->substr(0, pos);
    str->erase(0, pos);
    tree<SimGeneNode*>::leaf_iterator loc = tre->set_head(new SimGeneNode(root));
    while (!str->empty()) {
        char leftC = str->at(0); //left char must be one of the (, );
        size_t pos = str->find_first_of("()", 1);
        if (pos == std::string::npos) {
            loc = tre->parent(loc);
            break;
        }
        std::string locStr = str->substr(1, pos - 1);
        char rightC = str->at(pos); //right char must be one of the (, );
        str->erase(0, pos);
        if (leftC == '(') {
            auto locVec = splitStr(locStr);
            if (locVec.empty()) {
                error_exit("must have sth after1 '('!");
            }
            loc = tre->append_child(loc, new SimGeneNode(locVec.at(0)));
            for (int s = 1; s != locVec.size(); ++s) {
                loc = tre->insert_after(loc, new SimGeneNode(locVec.at(s)));
            }
        } else {
            loc = tre->parent(loc);
            auto locVec = splitStr(locStr);
            if (rightC == '(') {
                if (locVec.empty()) {
                    error_exit("must have loc between ')' and '('!");
                }
                for (auto & lo : locVec) {
                    loc = tre->insert_after(loc, new SimGeneNode(lo));
                }
            } else {
                for (auto & lo : locVec) {
                    loc = tre->insert_after(loc, new SimGeneNode(lo));
                }
            }
        }
    }
    return tre;
}

std::shared_ptr<tree<SimGeneNode*>> PhyloTree::buildTreePtrNode(std::string & db) {
    std::ifstream ifs(db.c_str());
    if(!ifs.is_open()) error_exit("can not open tree file: " + db);
    std::string line = "";
    std::queue<std::string> linQue;
    while(std::getline(ifs, line)){
        trimEnds(&line);
        if(line.empty()) continue;
        linQue.push(line);
    }
    ifs.close();
    if(mOptions->verbose) cerr << linQue.size() << " lines are readed for taxon tree" << "\n";
 
    std::shared_ptr<tree<SimGeneNode*>> tre(new tree<SimGeneNode*>());
    tre->set_head(new SimGeneNode("root"));
    auto locDummy = tre->append_child(tre->begin(), new SimGeneNode("dummy"));
    cCout("initiate a super root tree", tre->size());
    int numThreads = std::min(mOptions->thread, static_cast<int>(linQue.size()));
    std::thread consumerThreads[numThreads];
    for (int i = 0; i < numThreads; ++i) {
        consumerThreads[i] = std::thread([this, &i, &linQue, &tre, &locDummy]() {
            while (true) {
                std::unique_lock<std::mutex> lock(mtxTreR);
                if (linQue.empty()) {
                    lock.unlock();
                    break;
                }
                std::string linstr = linQue.front();
                linQue.pop();
                lock.unlock();
                if (linstr.empty() || linstr == "") continue;
                std::shared_ptr<tree<SimGeneNode*>> treTmp = buildTreeLoopPtrNode(&linstr);
                std::unique_lock<std::mutex> lock2(mtxTreW);
                locDummy = tre->insert_subtree(locDummy, treTmp->begin());
                lock2.unlock();
            }
        });
    }
    for (int i = 0; i < numThreads; ++i) {
        if (consumerThreads[i].joinable()) {
            consumerThreads[i].join();
        }
    }
    cCout("final tree depth", tre->max_depth(), 'g');
    return tre;
}

void PhyloTree::print_children_par(std::shared_ptr<tree<std::string *>>& tre, std::string fpath){
    std::ofstream *of = new std::ofstream();
    of->open(fpath.c_str(), std::ofstream::out);
    if(!of->is_open())
        error_exit("can not open " + fpath);
    *of << "par\tchildren\n";
    for (tree<std::string *>::post_order_iterator pos_loc = tre->begin_post(); pos_loc != tre->end_post(); pos_loc++){
        if(*(pos_loc.node->data) == "root" || *(pos_loc.node->data) == "dummy") continue;
        if(tre->is_valid(tre->parent(pos_loc))){
            *of << *(pos_loc.node->parent->data) << "\t" << *(pos_loc.node->data) << "\n";
        }
    }
    of->flush();
    of->close();
    if(of){
        delete of;
        of = nullptr;
    }
}
