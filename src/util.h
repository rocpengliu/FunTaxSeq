#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <sys/stat.h>
#include <algorithm>
#include <time.h>
#include <mutex>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <stdio.h>
#include <tuple>
#include <utility>
#include <numeric>
#include <string>  
#include <sstream>
#include <fstream>
#include <deque>
#include <regex>
#include <iomanip>
#include <cmath>
#include <cctype>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#include <limits.h>
#define GetCurrentDir getcwd
#endif

#ifdef WINDOWS
#include <windows.h>
#define GetExePath _getcwd
#else
#include <unistd.h>
#include <limits.h>
#define GetExePath getcwd
#endif

#include "options.h"

using namespace std;
extern mutex logmtx;

inline void loginfo(const string & s, bool next = true) {
    logmtx.lock();
    time_t tt = time(NULL);
    tm* t = localtime(&tt);
    if (next) {
        fprintf(stderr, "[\033[1;35m%02d:%02d:%02d\033[0m] %s\n", t->tm_hour, t->tm_min, t->tm_sec, s.c_str());
    } else {
        fprintf(stderr, "[\033[1;35m%02d:%02d:%02d\033[0m] %s\r", t->tm_hour, t->tm_min, t->tm_sec, s.c_str());
    }
    logmtx.unlock();
}

inline void loginfolong(const string & s) {
    logmtx.lock();
    time_t tt = time(NULL);
    tm* t = localtime(&tt);
    fprintf(stderr, "[%02d:%02d:%02d] %s\n", t->tm_hour, t->tm_min, t->tm_sec, s.c_str());
    logmtx.unlock();
}

template <typename T>
void ccCout(const T & str, bool first = true, bool last = true) {
    logmtx.lock();
    
    if(first){
        std::cout << str << "\t";
    } else {
        if (last) {
            std::cout << str << "\n";
        } else {
            std::cout << str << ";";
        }
    }
    
    logmtx.unlock();
}

template <typename T>
void cCout(const T & str, char color = 'r') {
    logmtx.lock();
    if (color == 'r') {
        std::cout << "\033[1;31m" << str << "\033[0m\n";
    } else if (color == 'g') {
        std::cout << "\033[1;32m" << str << "\033[0m\n";
    } else if (color == 'y') {
        std::cout << "\033[1;33m" << str << "\033[0m\n";
    } else if (color == 'b') {
        std::cout << "\033[1;34m" << str << "\033[0m\n";
    } else if (color == 'm') {
        std::cout << "\033[1;35m" << str << "\033[0m\n";
    } else if (color == 'c') {
        std::cout << "\033[1;36m" << str << "\033[0m\n";
    } else if (color == 'w') {
        std::cout << "\033[1;37m" << str << "\033[0m\n";
    } else {
        std::cout << "\033[1;30m" << str << "\033[0m\n";
    }
    logmtx.unlock();
}

template <typename T1, typename T2>
void cCout(const T1 & str1, const T2 & str2, char color = 'r') {
    logmtx.lock();
    if (color == 'r') {
        std::cout << "\033[1;31m" << str1 << " -> " << str2 << "\033[0m\n";
    } else if (color == 'g') {
        std::cout << "\033[1;32m" << str1 << " -> " << str2  << "\033[0m\n";
    } else if (color == 'y') {
        std::cout << "\033[1;33m" << str1 << " -> " << str2  << "\033[0m\n";
    } else if (color == 'b') {
        std::cout << "\033[1;34m" << str1 << " -> " << str2  << "\033[0m\n";
    } else if (color == 'm') {
        std::cout << "\033[1;35m" << str1 << " -> " << str2  << "\033[0m\n";
    } else if (color == 'c') {
        std::cout << "\033[1;36m" << str1 << " -> " << str2  << "\033[0m\n";
    } else if (color == 'w') {
        std::cout << "\033[1;37m" << str1 << " -> " << str2  << "\033[0m\n";
    } else {
        std::cout << "\033[1;30m" << str1 << " -> " << str2  << "\033[0m\n";
    }
    logmtx.unlock();
}

inline char complement(char base) {
    switch (base) {
        case 'A':
        case 'a':
            return 'T';
        case 'T':
        case 't':
            return 'A';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        default:
            return 'N';
    }
}

inline bool starts_with(string const & value, string const & start) {
    if (start.size() > value.size()) return false;
    return equal(start.begin(), start.end(), value.begin());
}

inline bool ends_with(string const & value, string const & ending) {
    if (ending.size() > value.size()) return false;
    return equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline string trim(const string& str) {
    string::size_type pos = str.find_first_not_of(' ');
    if (pos == string::npos) {
        return string("");
    }
    string::size_type pos2 = str.find_last_not_of(' ');
    if (pos2 != string::npos) {
        return str.substr(pos, pos2 - pos + 1);
    }
    return str.substr(pos);
}

inline void trimLeft(string &str, std::string sep = ";"){
    string::size_type pos = str.find_first_of(sep);
    if (pos != string::npos) {
        str.erase(0, pos + sep.length());
    }
}

inline void trimRight(string &str, std::string sep = ";"){
    if(str.back() == sep[0]){
        str.pop_back();
    }
}

inline void trimEnds(std::string* str) {
    if(str == nullptr) return;
    while (!str->empty() && (str->back() == '\n' || str->back() == '\r')) {
        str->pop_back();
    }
}

// Trim whitespace characters at the end of a char buffer
inline void trimEndsChar(char* buffer) {
    if (buffer == nullptr) return;
    // Find the length of the buffer
    size_t length = strlen(buffer);
    // Iterate from the end of the buffer backwards
    while (length > 0 && std::isspace(static_cast<unsigned char>(buffer[length - 1]))) {
        buffer[--length] = '\0'; // Replace whitespace with null terminator
    }
}

inline void removeEnds(std::string &str, char end = ';'){
    while (!str.empty() && str.back() == end) {
        str.pop_back();
    }
}

inline std::string removeExtension(std::string str, std::string extension){
    size_t pos = str.find(extension);
    if(pos == std::string::npos){
        return str;
    } else {
        return str.substr(0, pos);
    }
}
inline int splitStr(const string& str, vector<string> &ret_, string sep = "\t"){
    if (str.empty()) {
        return 0;
    }
    string tmp;
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;
    while (pos_begin != string::npos) {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos) {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        } else {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }
        ret_.push_back(tmp);
        tmp.clear();
    }
    ret_.shrink_to_fit();
    return 0;
}

inline std::vector<std::string> splitStr(string str, string sep = ";") {
    std::vector<std::string> ret_;
    if (str.empty()) {
        return ret_;
    }
    if(str.front() == sep[0]){
        str.erase(0, 1);
    }
    if(str.back() == sep[0]){
        str.pop_back();
    }
    int s = 0;
    for(const char & c : str){
        if(c == sep[0]){
            ++s;
        }
    }
    ret_.reserve(s + 1);
    string tmp = "";
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;
    while (pos_begin != string::npos) {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos) {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        } else {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }
        ret_.push_back(tmp);
        tmp.clear();
    }
    return ret_;
}

inline std::unordered_set<std::string> splitStr2(string str, string sep = ";") {
    std::unordered_set<std::string> ret_;
    if (str.empty()) {
        return ret_;
    }
    if(str.back() == sep[0]) {
        str.pop_back();
    }
    int s = 0;
    for(const char & c : str){
        if(c == sep[0]){
            ++s;
        }
    }
    ret_.reserve(s + 1);
    string tmp = "";
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;
    while (pos_begin != string::npos) {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos) {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        } else {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }
        ret_.insert(tmp);
        tmp.clear();
    }
    return ret_;
}

template<typename T>
inline std::string getStrVec(std::vector<T> vec, char sep = '|'){
    std::stringstream ss;
    for(const auto & it : vec){
        ss << it;
        if (&it != &(*vec.rbegin())) {
            ss << sep;
        }
    }
    return ss.str();
}

inline std::string getFirstNsSeps(std::string str, int n, char sep = ';'){
    std::size_t pos = str.find_first_of(sep);
    int i = 1;
    while(pos != std::string::npos){
        if(i == n){
            break;
        } else {
            pos = str.find_first_of(sep, pos + 1);
             ++i;
        }
    }
    if(i < n){
        return "";
    } else {
        return str.substr(0, pos);
    }
}
template <template<typename, typename...> class Container, typename T, typename... Args>
inline Container<T, Args...> splitStrInt(string str, std::string sep = ";") {
    Container<T, Args...> ret_;
    if (str.empty()) {
        return ret_;
    }
    if(str.front() == sep[0]){
        str.erase(0,1);
    }
    if(str.back() == sep[0]) {
        str.pop_back();
    }
    string tmp = "";
    string::size_type pos_begin = str.find_first_not_of(sep);
    string::size_type comma_pos = 0;
    while (pos_begin != string::npos) {
        comma_pos = str.find(sep, pos_begin);
        if (comma_pos != string::npos) {
            tmp = str.substr(pos_begin, comma_pos - pos_begin);
            pos_begin = comma_pos + sep.length();
        } else {
            tmp = str.substr(pos_begin);
            pos_begin = comma_pos;
        }
        //ret_.insert(static_cast<T>(std::stoul(tmp)));
        ret_.insert(tmp);
        tmp.clear();
    }
    return ret_;
}

inline int countFreq(std::string & s, std::string p){
    int M = s.length();
    int N = p.length();
    int res = 0;
    for(int i = 0; i <= M - N; i ++){
        int j;
        for(j = 0; j < N; j++){
            if(s[i+j] != p[j]) break;
        }
        if(j == N){
            ++res;
            j = 0;
        }
    }
    return res;
}

inline string replace(const string& str, const string& src, const string& dest) {
    string ret;
    string::size_type pos_begin = 0;
    string::size_type pos = str.find(src);
    while (pos != string::npos) {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        ret += dest;
        pos_begin = pos + 1;
        pos = str.find(src, pos_begin);
    }
    if (pos_begin < str.length()) {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline string reverse(const string& str) {
    string ret(str.length(), 0);
    for (int pos = 0; pos < str.length(); pos++) {
        ret[pos] = str[str.length() - pos - 1];
    }
    return ret;
}

inline string basename(string filename) {
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos)
        return filename;
    else if (pos == filename.length() - 1)
        return ""; // a bad filename
    else
        return filename.substr(pos + 1, filename.length() - pos - 1);
}

inline string get_current_dir() {
    char buff[FILENAME_MAX]; //create string buffer to hold path
    GetCurrentDir(buff, FILENAME_MAX);
    string current_working_dir(buff);
    return current_working_dir;
}

inline string get_upper_dir() {
    std::string cwd = get_current_dir();
    std::string::size_type bepos = cwd.find_last_of("/");
    cwd.erase(bepos);
    return (cwd);
}

#ifdef WINDOWS

std::string getexepath() {
    char result[ MAX_PATH ];
    return std::string(result, GetModuleFileName(NULL, result, MAX_PATH));
}
#else

inline string GetExePath() {
    char result[ PATH_MAX ];
    ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
    return std::string(result, (count > 0) ? count : 0);
}

inline string get_funtaxseq_dir() {
    std::string cwd = GetExePath();
    std::string::size_type bepos = cwd.find("/bin");
    cwd.erase(bepos);
    return (cwd);
}
#endif

inline string dirname(const string& filename) {
    string::size_type pos = filename.find_last_of('/');
    if (pos == string::npos) {
        return "./";
    } else {
        return filename.substr(0, pos + 1);
    }
}

inline string joinpath(const string& dirname, const string& basename) {
    if (dirname[dirname.length() - 1] == '/') {
        return dirname + basename;
    } else {
        return dirname + "/" + basename;
    }
}

inline string checkDirEnd(const string & dirname) {
    return (dirname[dirname.length() - 1] == '/' ? dirname : dirname + '/');
}

//Check if a string is a file or directory

inline bool file_exists(const string& s) {
    bool exists = false;
    if (s.length() > 0) {
        struct stat status;
        int result = stat(s.c_str(), &status);
        if (result == 0) {
            exists = true;
        }
    }
    return exists;
}

// check if a string is a directory

inline bool is_directory(const string& path) {
    bool isdir = false;
    struct stat status;
    // visual studion use _S_IFDIR instead of S_IFDIR
    // http://msdn.microsoft.com/en-us/library/14h5k7ff.aspx
#ifdef _MSC_VER
#define S_IFDIR _S_IFDIR
#endif
    stat(path.c_str(), &status);
    if (status.st_mode & S_IFDIR) {
        isdir = true;
    }
    // #endif
    return isdir;
}

inline void check_file_valid(const string& s) {
    if (!file_exists(s)) {
        cerr << "ERROR: file '" << s << "' doesn't exist, quit now" << endl;
        exit(-1);
    }
    if (is_directory(s)) {
        cerr << "ERROR: '" << s << "' is a folder, not a file, quit now" << endl;
        exit(-1);
    }
}

inline bool check_filename_valid(const string& s) {
    return 0 < trim(s).length() && trim(s).length() <= 255 && regex_match(s, regex("^[A-Za-z0-9_\\.\\-]+$"));
}

inline void check_file_writable(const string& s) {
    string dir = dirname(s);
    if (!file_exists(dir)) {
        cerr << "ERROR: '" << dir << " doesn't exist. Create this folder and run this command again." << endl;
        exit(-1);
    }
    if (is_directory(s)) {
        cerr << "ERROR: '" << s << "' is not a writable file, quit now" << endl;
        exit(-1);
    }
}

// Remove non alphabetic characters from a string

inline string str_keep_alpha(const string& s) {
    string new_str;
    for (size_t it = 0; it < s.size(); it++) {
        if (isalpha(s[it])) {
            new_str += s[it];
        }
    }
    return new_str;
}

// Remove invalid sequence characters from a string
inline void str_keep_valid_sequence(string& s, bool forceUpperCase = false) {
    size_t total = 0;
    const char case_gap = 'a' - 'A';
    for (size_t it = 0; it < s.size(); it++) {
        char c = s[it];
        if (forceUpperCase && c >= 'a' && c <= 'z') {
            c -= case_gap;
        }
        if (isalpha(c) || c == '-' || c == '*') {
            s[total] = c;
            ++total;
        }
    }

    s.resize(total);
}

inline int find_with_right_pos(const string& str, const string& pattern, int start = 0) {
    int pos = str.find(pattern, start);
    if (pos < 0)
        return -1;
    else
        return pos + pattern.length();
}

inline void str2upper(string& s) {
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))toupper);
}

inline void str2lower(string& s) {
    transform(s.begin(), s.end(), s.begin(), (int (*)(int))tolower);
}

inline char num2qual(int num) {
    if (num > 127 - 33)
        num = 127 - 33;
    if (num < 0)
        num = 0;

    char c = num + 33;
    return c;
}

inline void error_exit(const string& msg) {
    cerr << "ERROR: " << msg << endl;
    exit(-1);
}

inline string trimName(string& str) {
    // string strnew;
    str.erase(str.begin());
    string suffixStartCharacters = " /\t\r";
    size_t n = str.find_first_of(suffixStartCharacters);
    if (n != string::npos) {
        return str.erase(n);
    } else {
        return str;
    }
}

inline string trimName2(string src) {
    // string strnew;
    std::string str = src;
    string suffixStartCharacters = " /\t\r";
    size_t n = str.find_first_of(suffixStartCharacters);
    if (n != string::npos) {
        return str.erase(n);
    } else {
        return str;
    }
}

inline string removeStr(const string &str, const string &src) {
    string ret;
    string::size_type pos_begin = 0;
    string::size_type pos = str.find(src);
    while (pos != string::npos) {
        ret.append(str.data() + pos_begin, pos - pos_begin);
        pos_begin = pos + src.length();
        pos = str.find(src, pos_begin);
    }
    if (pos_begin < str.length()) {
        ret.append(str.begin() + pos_begin, str.end());
    }
    return ret;
}

inline string removeStrs(const string &str) {
    string ret;
    string::size_type pos_begin = 0;
    string::size_type pos = str.find("\tK");
    if (pos != string::npos) {
        ret.append(str.data(), pos);
    } else {
        pos = str.find("\tU");
        if (pos != string::npos) {
            ret.append(str.data(), pos);
        } else {
            return str;
        }
    }
    return ret;
}

template<typename T1, typename T2>
double getPer(T1 v1, T2 v2, bool per = true) {
    double ret = 0.00;
    if(v1 == 0)
        return ret;
    if(per){
        ret = std::round((static_cast<double>(v1 * 100.0)/v2) * 10000.0) / 10000.0;
    } else {
        ret = std::round((static_cast<double>(v1)/v2) * 10000.0) / 10000.0;
    }
    return ret;
}

template<typename T>
std::string unkown2Str(const T & t) {
    std::ostringstream os;
    os << t;
    return os.str();
}

inline void getUniqVec(std::vector<std::string> &orgVec) {
    std::unordered_set<std::string> s;
    std::vector<std::string>::iterator itr = orgVec.begin();
    for (auto curr = orgVec.begin(); curr != orgVec.end(); ++curr) {
        if (s.insert(*curr).second)
            *itr++ = *curr;
    }
    orgVec.erase(itr, orgVec.end());
}

inline void stripChar(std::string & s) {
    for (auto it = s.begin(); it != s.end(); ++it) {
        if (!isalpha(*it)) {
            s.erase(it);
            it--;
        }
    }
}

template<typename T>
inline std::string convertSeconds(const T & t) {
    long seconds, minutes, hours;
    seconds = long(t);
    minutes = t / 60;
    hours = minutes / 60;
    std::stringstream ss;
    ss << hours << " hours " << minutes % 60 << " minutes " << seconds % 60 << " seconds";
    return ss.str();
}

template<typename T>
inline T getMapMaxKey(std::map<T, int>&m){
    int maxValKey = 0;
    T key;
    for(const auto & it : m){
        if(it.second > maxValKey){
            key = it.first;
            maxValKey = it.second;
        }
    }
    return key;
}
#endif /* UTIL_H */
