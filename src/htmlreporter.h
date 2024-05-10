#ifndef HTML_REPORTER_H
#define HTML_REPORTER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <utility>
#include <set>
#include <algorithm>
#include "options.h"
#include "stats.h"
#include "filterresult.h"

using namespace std;

class HtmlReporter{
public:
    HtmlReporter(Options* opt);
    ~HtmlReporter();
    void setDupHist(int* dupHist, double* dupMeanGC, double dupRate);
    void setInsertHist(long* insertHist, int insertSizePeak);
    void report(FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2 = NULL, Stats* postStats2 = NULL);

    static void outputRow(ofstream& ofs, string key, long value);
    static void outputRow(ofstream& ofs, string key, string value);
    static string formatNumber(long number);
    static string getPercents(long numerator, long denominator);
    static std::string highligher(std::string & str, std::map<int, std::string> & snpsMap);
    static std::string highligher(std::string & str, std::set<int> & snpsSet);
   
    
private:
    const string getCurrentSystemTime();
    void printHeader(ofstream& ofs);
    void printCSS(ofstream& ofs);
    void printJS(ofstream& ofs);
    void printFooter(ofstream& ofs);
    void reportDuplication(ofstream& ofs);
    void reportInsertSize(ofstream& ofs, int isizeLimit);
    void printSummary(ofstream& ofs, FilterResult* result, Stats* preStats1, Stats* postStats1, Stats* preStats2, Stats* postStats2);
    
private:
    Options* mOptions;
    int* mDupHist;
    double* mDupMeanGC;
    double mDupRate;
    long* mInsertHist;
    int mInsertSizePeak;
};


#endif