#include "phylotree.h"

PhyloTree::PhyloTree(PhyloOptions *& mOptions) {
    this->mOptions = mOptions;
    geneTree = nullptr;
    geneNodeTree = nullptr;
    taxonTree = nullptr;
    geneAnoMap.clear();
    orthAnoMap.clear();
    geneDNADupMap.clear();
    geneProDupMap.clear();
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

    std::thread dFTreThread = std::thread([this](){
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
    std::thread readGeno = std::thread([this, &geneAnnoQueue]{
        geneAnnoQueue = readGZ(mOptions->geneAno);
    });

    std::queue<std::string> orthAnnoQueue;
    std::thread readOrth = std::thread([this, &orthAnnoQueue]{
       orthAnnoQueue = readGZ(mOptions->orthAno);
    });

    std::queue<std::string> geneDNADupQueue;
    std::thread geneDNADup = std::thread([this, &geneDNADupQueue]{
        geneDNADupQueue = readGZ(mOptions->geneDNADup);
    });

    std::queue<std::string> geneProDupQueue;
    std::thread geneProDup = std::thread([this, &geneProDupQueue]{
        geneProDupQueue = readGZ(mOptions->geneProDup);
    });

    if(readGeno.joinable()){
        readGeno.join();
    }
    if(readOrth.joinable()){
        readOrth.join();
    }
    if(geneDNADup.joinable()){
        geneDNADup.join();
    }
    if(geneProDup.joinable()){
        geneProDup.join();
    }
    readGeneAnno(geneAnnoQueue);
    readOrthAnno(orthAnnoQueue);
    readGeneDup(geneDNADupQueue, 'd');
    readGeneDup(geneProDupQueue, 'p');
    if (mOptions->verbose)
        loginfo("start to build taxon tree!");
    taxonTree = buildTreePtr(mOptions->tTree, "taxon");
    if(mOptions->verbose) loginfo("finished to build taxon tree!");
    if(taxonTree->size() < 3) error_exit("built taxon tree size must be no less than 2: ");
    if(mOptions->verbose) cerr << "taxon tree size is " << taxonTree->size() << " and has " << taxonTree->begin().number_of_descent() << " descents" << "\n";
    
    if(mOptions->verbose) loginfo("start to build gene ortholog tree!");
    //geneNodeTree = buildTreePtrNode(mOptions->gTree);
    //populateGeneTre();
    geneTree = buildTreePtr(mOptions->gTree, "ortholog");
    geneHeadIterator = std::find_if(geneTree->begin(),)
    if(mOptions->verbose) loginfo("finished to build gene ortholog tree!");
    if(geneTree->size() < 1) error_exit("built gene tree size must be no less than 1: ");
    if(mOptions->verbose) cerr << "gene ortholog tree size is " << geneTree->size() << " and has " << geneTree->begin().number_of_descent() << " descents" << "\n";

//     tree<std::string*>::post_order_iterator locf;
//     locf = std::find_if(geneTree->begin_post(), geneTree->end_post(),
//                 [](std::string* & itp) {
//                     return *itp == "9807975at2";
//                 });
//     cCout("searching in a tree");
//     if (geneTree->is_valid(locf)) {
//         cCout("locf leaf is", *(locf.node->data), 'y');
//     }
//     for (tree<std::string*>::pre_order_iterator it = geneTree->begin(); it != geneTree->end(); it++) {
//        //if(tre.is_valid(it)){
//        if (it.node->data) {
//         cCout(*(it.node->data));
//        }
//    }
}

/*
void PhyloTree::init(){
    std::ifstream ifs(mOptions->tTree.c_str());
    if(!ifs.is_open()) error_exit("can not open taxon tree: " + mOptions->tTree);
    
    std::string line = "";
    std::queue<std::string> linQue;
    while(std::getline(ifs, line)){
        trimEnds(&line);
        if(line.empty()) continue;
        linQue.push(line);
    }
    ifs.close();
    if(mOptions->verbose) cerr << linQue.size() << " lines are readed for taxon tree" << "\n";
    
    if(mOptions->verbose) loginfo("start to build taxon tree!");
    taxonTree = buildTreePtr(linQue, mOptions->thread);
    //taxonTree = buildTreeIntPtr(linQue, mOptions->thread);
    if(mOptions->verbose) loginfo("finished to build taxon tree!");
    if(taxonTree->size() < 3) error_exit("built taxon tree size must be no less than 2: ");
    if(mOptions->verbose) cerr << "taxon tree size is " << taxonTree->size() << " and has " << taxonTree->begin().number_of_descent() << " descents" << "\n";
    ifs.open(mOptions->gTree.c_str());
    if(!ifs.is_open()) error_exit("can not open gene ortholog tree: " + mOptions->gTree);
    //std::cin.get();

    line = "";
    while(!linQue.empty()) linQue.pop();
    while(std::getline(ifs, line)){
        trimEnds(&line);
        if(line.empty()) continue;
        linQue.push(line);
    }
    ifs.close();
    if(mOptions->verbose) cerr << linQue.size() << " lines are readed for gene ortholog tree" << "\n";
    
    if(mOptions->verbose) loginfo("start to build gene ortholog tree!");
    geneTree = buildTreePtr(linQue, mOptions->thread);
    if(mOptions->verbose) loginfo("finished to build gene ortholog tree!");
    if(geneTree->size() < 1) error_exit("built gene tree size must be no less than 1: ");
    if(mOptions->verbose) cerr << "gene ortholog tree size is " << geneTree->size() << " and has " << geneTree->begin().number_of_descent() << " descents" << "\n";
    //std::cin.get();

    ifs.open(mOptions->othMap.c_str());
    if(!ifs.is_open()) error_exit("can not open gene ortholog tree: " + mOptions->othMap);
    line = "";
    std::vector<std::string> strVec;
    while(std::getline(ifs, line)){
        trimEnds(&line);
        splitStr(line, strVec);
        if(strVec.size() == 2){
           orthMap[strVec[1]] = strVec[0];
        }
        strVec.clear();
        if(mOptions->othMap.size() % 5000000 == 0) {
            if(mOptions->verbose) cerr << "reading ogmap  " << orthMap.size() << "\n";
        }
    }
    ifs.close();
    if(mOptions->verbose) cerr << "gene ortholog map is " << orthMap.size() << "\n";
    //std::cin.get();
}

*/

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
void PhyloTree::populateGeneTre(){
    if(mOptions->verbose) loginfo("populate gene tree");
    tree<SimGeneNode *>::iterator it;
    for (it = geneNodeTree->begin(); it != geneNodeTree->end(); ++it){
        auto locn = it.node->data->id;
        auto loc = orthAnoMap.find(locn);
        if(loc == orthAnoMap.end())
            continue;
        it.node->data->anno = loc->second->anno;
        it.node->data->goSet = loc->second->goSet;
        it.node->data->koSet = loc->second->koSet;;
    }
    if(mOptions->verbose) loginfo("populate gene tree done");

    //delete the orthmap to save memeory
}

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

void PhyloTree::readGeneAnno(std::queue<std::string>& geneAnnoQueue){
    if(mOptions->verbose){
        loginfo("start to populate gene map");
    }
    int numThreads = std::min(static_cast<int>(geneAnnoQueue.size()), mOptions->thread);
    std::thread consumerThreads[numThreads];
    for(int i = 0; i < numThreads; ++i){
        consumerThreads[i] = std::thread([this, &geneAnnoQueue, &i](){
            while(true){
                std::unique_lock<std::mutex> lock(mtxTreR);
                if(geneAnnoQueue.empty()){
                    lock.unlock();
                    break;
                }
                std::string line = geneAnnoQueue.front();
                geneAnnoQueue.pop();
                if(geneAnnoQueue.size() >= 100000 && geneAnnoQueue.size() % 100000 == 0){
                    std::string msg = "gene anno remains: " + std::to_string(geneAnnoQueue.size());
                    loginfo(msg, false);
                }
                lock.unlock();
                if(line.empty())
                    continue;
                std::vector<std::string> strVec;
                splitStr(line, strVec);
                if(strVec.size() == 6){
                    GeneNode *tmp = new GeneNode();
                    tmp->id = strVec[0];
                    tmp->par = strVec[1];
                    tmp->taxon = strVec[2];
                    tmp->anno = strVec[3];
                    if(strVec[4] != "0"){
                        tmp->goSet = splitStrInt<std::set, std::string>(strVec[4]);
                    }
                    if (strVec[5] != "0"){
                         tmp->koSet = splitStrInt<std::set, std::string>(strVec[5]);
                    }
                    std::unique_lock<std::mutex> lock2(mtxTreW);
                    geneAnoMap[strVec[0]] = tmp;
                    lock2.unlock();
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
        loginfo("read gene map done with size " + std::to_string(geneAnoMap.size()));
    }
}

void PhyloTree::readOrthAnno(std::queue<std::string>& orthAnnoQueue){
    if(mOptions->verbose){
        loginfo("start to populate ortho map");
    }
    int numThreads = std::min(static_cast<int>(orthAnnoQueue.size()), mOptions->thread);
    std::thread consumerThreads[numThreads];
    for(int i = 0; i < numThreads; ++i){
        consumerThreads[i] = std::thread([this, &orthAnnoQueue, &i](){
            while(true){
                std::unique_lock<std::mutex> lock(mtxTreR);
                if(orthAnnoQueue.empty()){
                    lock.unlock();
                    break;
                }
                std::string line = orthAnnoQueue.front();
                orthAnnoQueue.pop();
                if(orthAnnoQueue.size() >= 100000 && orthAnnoQueue.size() % 100000 == 0){
                    std::string msg = "ortholog remains: " + std::to_string(orthAnnoQueue.size());
                    loginfo(msg, false);
                }
                lock.unlock();
                if(line.empty())
                    continue;
                std::vector<std::string> strVec;
                splitStr(line, strVec);
                if(strVec.size() == 6){
                    GeneNode *tmp = new GeneNode();
                    tmp->id = strVec[0];
                    tmp->par = strVec[1];
                    tmp->taxon = strVec[2];
                    tmp->anno = strVec[3];
                    if(strVec[4] != "0"){
                        tmp->goSet = splitStrInt<std::set, std::string>(strVec[4], ";");
                    }
                    if (strVec[5] != "0"){
                        tmp->koSet = splitStrInt<std::set, std::string>(strVec[5], ";");
                    }
                    std::unique_lock<std::mutex> lock2(mtxTreW);
                    orthAnoMap[strVec[0]] = tmp;
                    lock2.unlock();
                }
            } });
    }
    for(int i = 0; i < numThreads; i++){
        if(consumerThreads[i].joinable()){
            consumerThreads[i].join();
        }
    }
    if(mOptions->verbose){
        loginfo("populate ortho map done with size " + std::to_string(orthAnoMap.size()));
    }
}

void PhyloTree::readGeneDup(std::queue<std::string>& geneDupQueue, char type){
    if(mOptions->verbose){
        loginfo("start to populate gene dup map");
    }
    int numThreads = std::min(static_cast<int>(geneDupQueue.size()), mOptions->thread);
    std::thread consumerThreads[numThreads];
    for(int i = 0; i < numThreads; ++i){
        consumerThreads[i] = std::thread([this, &geneDupQueue, &i, &type](){
            while(true){
                std::unique_lock<std::mutex> lock(mtxTreR);
                if(geneDupQueue.empty()){
                    lock.unlock();
                    break;
                }
                std::string line = geneDupQueue.front();
                geneDupQueue.pop();
                if(geneDupQueue.size() >= 100000 && geneDupQueue.size() % 100000 == 0){
                    std::string msg = "gene dup remains: " + std::to_string(geneDupQueue.size());
                    loginfo(msg, false);
                }
                lock.unlock();
                if(line.empty())
                    continue;
                std::vector<std::string> strVec;
                splitStr(line, strVec);
                if(strVec.size() == 2){
                    std::unique_lock<std::mutex> lock2(mtxTreW);
                    if(type == 'd'){
                        geneDNADupMap[strVec.at(0)] = strVec.at(1);
                    } else if(type == 'p') {
                        geneProDupMap[strVec.at(0)] = strVec.at(1);
                    }
                    lock2.unlock();
                }
            } });
    }
    for(int i = 0; i < numThreads; i++){
        if(consumerThreads[i].joinable()){
            consumerThreads[i].join();
        }
    }
    if(mOptions->verbose){
        loginfo("populate dup map done with size " + std::to_string(geneDNADupMap.size()));
    }
}

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
