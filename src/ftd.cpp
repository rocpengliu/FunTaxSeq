#include <stdio.h>
#include <time.h>
#include <sstream>
#include <string>

#include "util.h"
#include "unittest.h"
#include "cmdline.h"
#include "funtaxoptions.h"
#include "funtaxdecoder.h"

using namespace std;

string command;
mutex logmtx;
std::mutex mtxTreR;
std::mutex mtxTreW;

int main(int argc, char* argv[]) {
    // display version info if no argument is given
    if (argc == 1) {
        cerr << "FunTaxDecoder: decoding the output of FunTaxSeq to taxonic and functional profilings." << endl << "version " << FUNTAXSEQ_VER << endl;
    }
    if (argc == 2 && (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0)) {
        cerr << "funtaxseq " << FUNTAXSEQ_VER << endl;
        return 0;
    }
    cmdline::parser cmd;
    cmd.add<string>("samdir", 's', "sample directory", false, "");
    cmd.add<string>("outprefix", 'o', "output file prefix", false, "");
    cmd.add<string>("otre", 't', "output tree par kid file", false, "");
    cmd.add<string>("database", 'b', "datbase directory", false, "");
    cmd.add("useogtree", 'g', "If specified, using ortholog tree(slow)");
    cmd.add<int>("thread", 'w', "worker thread number, default is 4", false, 4);
    cmd.add("debug", 'd', "If specified, print debug");
    cmd.add("verbose", 'V', "output verbose");

    cmd.parse_check(argc, argv);
    if (argc == 1) {
        cerr << cmd.usage() << endl;
        return 0;
    }
    time_t t11 = time(NULL);
    stringstream ss;
    for (int i = 0; i < argc; i++) {
        ss << argv[i] << " ";
    }

    command = ss.str();
    PhyloOptions * opt = new PhyloOptions();
    opt->sampleDir = cmd.get<string>("samdir");
    opt->db = cmd.get<string>("database");
    opt->prefix = cmd.get<string>("outprefix");
    opt->outTaxon =  opt->prefix + "_taxon_abundance.txt";
    opt->outFun = opt->prefix + "_raw_func_abundance.txt";
    opt->outPureFun = opt->prefix + "_gene_go_ko_func_abundance.txt";
    opt->outGeneFun = opt->prefix + "_gene_func_abundance.txt";
    if(cmd.exist("useogtree")){
        opt->gTree = joinpath(opt->db, "ogs_tree.tre");
    }
    opt->tTree = joinpath(opt->db, "taxon_tree.tre");
    opt->taxonGenomeSize = joinpath(opt->db, "taxon_genome_size.tab.gz");
    //opt->orthGeneSize = joinpath(opt->db, "orth_gene_size.tab.gz");
    opt->geneAno = joinpath(opt->db, "genes_anno.tab.gz");
    opt->orthAno = joinpath(opt->db, "ogs_anno.tab.gz");
    opt->geneDNADup = joinpath(opt->db, "dna_duplicated.tab.gz");
    opt->geneProDup = joinpath(opt->db, "pro_duplicated.tab.gz");
    opt->verbose = cmd.exist("verbose");
    opt->debug = cmd.exist("debug");
    opt->thread = cmd.get<int>("thread");
    opt->outTree = cmd.get<string>("otre");
    opt->parseSample();

    FunTaxDecoder * ftDecoder = new FunTaxDecoder(opt);
    ftDecoder->process();
    if(ftDecoder){
        delete ftDecoder;
        ftDecoder = NULL;
    }
    time_t t22 = time(NULL);
    cerr << endl << command << endl;
    cerr << "funtaxdecoder v" << FUNTAXSEQ_VER << " time used: " << convertSeconds((t22) - t11) << " to process " << opt->samples.size() << " samples" << endl;
    if(opt){
        delete opt;
        opt = nullptr;
    }
    return 0;
}