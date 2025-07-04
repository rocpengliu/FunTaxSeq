#include "transsearcher.h"

TransSearcher::TransSearcher(Options * & opt, BwtFmiDB* & mBwtfmiDB) {
    mOptions = opt;
    tbwtfmiDB = mBwtfmiDB;
    min_match_length = 0;
    //match_str.clear();
    blosum_subst = {
        {'A',
         {'S', 'V', 'T', 'G', 'C', 'P', 'M', 'K', 'L', 'I', 'E', 'Q', 'R', 'Y', 'F', 'H', 'D', 'N', 'W'}},
        {'R',
         {'K', 'Q', 'H', 'E', 'N', 'T', 'S', 'M', 'A', 'Y', 'P', 'L', 'G', 'D', 'V', 'W', 'F', 'I', 'C'}},
        {'N',
         {'S', 'H', 'D', 'T', 'K', 'G', 'E', 'Q', 'R', 'Y', 'P', 'M', 'A', 'V', 'F', 'L', 'I', 'C', 'W'}},
        {'D',
         {'E', 'N', 'S', 'Q', 'T', 'P', 'K', 'H', 'G', 'R', 'A', 'V', 'Y', 'F', 'M', 'I', 'C', 'W', 'L'}},
        {'C',
         {'A', 'V', 'T', 'S', 'M', 'L', 'I', 'Y', 'W', 'F', 'P', 'K', 'H', 'G', 'Q', 'D', 'N', 'R', 'E'}},
        {'Q',
         {'E', 'K', 'R', 'S', 'M', 'H', 'D', 'N', 'Y', 'T', 'P', 'A', 'V', 'W', 'L', 'G', 'F', 'I', 'C'}},
        {'E',
         {'Q', 'D', 'K', 'S', 'H', 'N', 'R', 'T', 'P', 'A', 'V', 'Y', 'M', 'G', 'W', 'F', 'L', 'I', 'C'}},
        {'G',
         {'S', 'N', 'A', 'D', 'W', 'T', 'P', 'K', 'H', 'E', 'Q', 'R', 'V', 'Y', 'F', 'M', 'C', 'L', 'I'}},
        {'H',
         {'Y', 'N', 'E', 'Q', 'R', 'S', 'F', 'K', 'D', 'W', 'T', 'P', 'M', 'G', 'A', 'V', 'L', 'I', 'C'}},
        {'I',
         {'V', 'L', 'M', 'F', 'Y', 'T', 'C', 'A', 'S', 'W', 'P', 'K', 'H', 'E', 'Q', 'D', 'N', 'R', 'G'}},
        {'L',
         {'M', 'I', 'V', 'F', 'Y', 'T', 'C', 'A', 'W', 'S', 'K', 'Q', 'R', 'P', 'H', 'E', 'N', 'G', 'D'}},
        {'K',
         {'R', 'E', 'Q', 'S', 'N', 'T', 'P', 'M', 'H', 'D', 'A', 'V', 'Y', 'L', 'G', 'W', 'F', 'I', 'C'}},
        {'M',
         {'L', 'V', 'I', 'F', 'Q', 'Y', 'W', 'T', 'S', 'K', 'C', 'R', 'A', 'P', 'H', 'E', 'N', 'G', 'D'}},
        {'F',
         {'Y', 'W', 'M', 'L', 'I', 'V', 'H', 'T', 'S', 'C', 'A', 'K', 'G', 'E', 'Q', 'D', 'N', 'R', 'P'}},
        {'P',
         {'T', 'S', 'K', 'E', 'Q', 'D', 'A', 'V', 'M', 'H', 'G', 'N', 'R', 'Y', 'L', 'I', 'C', 'W', 'F'}},
        {'S',
         {'T', 'N', 'A', 'K', 'G', 'E', 'Q', 'D', 'P', 'M', 'H', 'C', 'R', 'V', 'Y', 'F', 'L', 'I', 'W'}},
        {'T',
         {'S', 'V', 'N', 'A', 'P', 'M', 'K', 'L', 'I', 'E', 'Q', 'C', 'D', 'R', 'Y', 'W', 'F', 'H', 'G'}},
        {'W',
         {'Y', 'F', 'M', 'T', 'L', 'H', 'G', 'Q', 'C', 'V', 'S', 'K', 'I', 'E', 'R', 'A', 'P', 'D', 'N'}},
        {'Y',
         {'F', 'W', 'H', 'V', 'M', 'L', 'I', 'Q', 'T', 'S', 'K', 'E', 'C', 'N', 'R', 'A', 'P', 'G', 'D'}},
        {'V',
         {'I', 'M', 'L', 'T', 'A', 'Y', 'F', 'C', 'S', 'P', 'K', 'E', 'Q', 'W', 'H', 'G', 'D', 'N', 'R'}}};

    std::memset(nuc2int, std::numeric_limits<uint8_t>::max(), sizeof (nuc2int));
    nuc2int['A'] = nuc2int['a'] = 0;
    nuc2int['C'] = nuc2int['c'] = 1;
    nuc2int['G'] = nuc2int['g'] = 2;
    nuc2int['T'] = nuc2int['t'] = 3;
    nuc2int['U'] = nuc2int['u'] = 3;
    std::memset(compnuc2int, std::numeric_limits<uint8_t>::max(), sizeof (compnuc2int));
    compnuc2int['A'] = compnuc2int['a'] = 3;
    compnuc2int['C'] = compnuc2int['c'] = 2;
    compnuc2int['G'] = compnuc2int['g'] = 1;
    compnuc2int['T'] = compnuc2int['t'] = 0;
    compnuc2int['U'] = compnuc2int['u'] = 0;

    std::memset(aa2int, std::numeric_limits<uint8_t>::min(), sizeof (aa2int));
    aa2int['A'] = 0;
    aa2int['R'] = 1;
    aa2int['N'] = 2;
    aa2int['D'] = 3;
    aa2int['C'] = 4;
    aa2int['Q'] = 5;
    aa2int['E'] = 6;
    aa2int['G'] = 7;
    aa2int['H'] = 8;
    aa2int['I'] = 9;
    aa2int['L'] = 10;
    aa2int['K'] = 11;
    aa2int['M'] = 12;
    aa2int['F'] = 13;
    aa2int['P'] = 14;
    aa2int['S'] = 15;
    aa2int['T'] = 16;
    aa2int['W'] = 17;
    aa2int['Y'] = 18;
    aa2int['V'] = 19;
    blosum62diag[aa2int['A']] = 4;
    blosum62diag[aa2int['R']] = 5;
    blosum62diag[aa2int['N']] = 6;
    blosum62diag[aa2int['D']] = 6;
    blosum62diag[aa2int['C']] = 9;
    blosum62diag[aa2int['Q']] = 5;
    blosum62diag[aa2int['E']] = 5;
    blosum62diag[aa2int['G']] = 6;
    blosum62diag[aa2int['H']] = 8;
    blosum62diag[aa2int['I']] = 4;
    blosum62diag[aa2int['L']] = 4;
    blosum62diag[aa2int['K']] = 5;
    blosum62diag[aa2int['M']] = 5;
    blosum62diag[aa2int['F']] = 6;
    blosum62diag[aa2int['P']] = 7;
    blosum62diag[aa2int['S']] = 4;
    blosum62diag[aa2int['T']] = 5;
    blosum62diag[aa2int['W']] = 11;
    blosum62diag[aa2int['Y']] = 7;
    blosum62diag[aa2int['V']] = 4;

    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'R']] = -1;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'N']] = -2;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'D']] = -2;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'C']] = 0;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'Q']] = -1;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'E']] = -1;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'G']] = 0;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'H']] = -2;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'I']] = -1;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'L']] = -1;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'K']] = -1;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'M']] = -1;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'F']] = -2;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'P']] = -1;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'S']] = 1;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'T']] = 0;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'W']] = -3;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'Y']] = -2;
    b62[aa2int[(uint8_t) 'A']][aa2int[(uint8_t) 'V']] = 0;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'A']] = -1;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'N']] = 0;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'D']] = -2;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'C']] = -3;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'Q']] = 1;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'E']] = 0;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'G']] = -2;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'H']] = 0;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'I']] = -3;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'L']] = -2;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'K']] = 2;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'M']] = -1;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'F']] = -3;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'P']] = -2;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'S']] = -1;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'W']] = -3;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'Y']] = -2;
    b62[aa2int[(uint8_t) 'R']][aa2int[(uint8_t) 'V']] = -3;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'A']] = -2;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'R']] = 0;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'D']] = 1;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'C']] = -3;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'Q']] = 0;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'E']] = 0;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'G']] = 0;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'H']] = 1;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'I']] = -3;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'L']] = -3;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'K']] = 0;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'M']] = -2;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'F']] = -3;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'P']] = -2;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'S']] = 1;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'T']] = 0;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'W']] = -4;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'Y']] = -2;
    b62[aa2int[(uint8_t) 'N']][aa2int[(uint8_t) 'V']] = -3;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'A']] = -2;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'R']] = -2;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'N']] = 1;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'C']] = -3;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'Q']] = 0;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'E']] = 2;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'G']] = -1;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'H']] = -1;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'I']] = -3;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'L']] = -4;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'K']] = -1;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'M']] = -3;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'F']] = -3;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'P']] = -1;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'S']] = 0;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'W']] = -4;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'Y']] = -3;
    b62[aa2int[(uint8_t) 'D']][aa2int[(uint8_t) 'V']] = -3;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'A']] = 0;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'R']] = -3;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'N']] = -3;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'D']] = -3;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'Q']] = -3;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'E']] = -4;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'G']] = -3;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'H']] = -3;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'I']] = -1;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'L']] = -1;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'K']] = -3;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'M']] = -1;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'F']] = -2;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'P']] = -3;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'S']] = -1;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'W']] = -2;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'Y']] = -2;
    b62[aa2int[(uint8_t) 'C']][aa2int[(uint8_t) 'V']] = -1;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'A']] = -1;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'R']] = 1;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'N']] = 0;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'D']] = 0;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'C']] = -3;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'E']] = 2;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'G']] = -2;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'H']] = 0;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'I']] = -3;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'L']] = -2;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'K']] = 1;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'M']] = 0;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'F']] = -3;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'P']] = -1;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'S']] = 0;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'W']] = -2;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'Y']] = -1;
    b62[aa2int[(uint8_t) 'Q']][aa2int[(uint8_t) 'V']] = -2;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'A']] = -1;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'R']] = 0;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'N']] = 0;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'D']] = 2;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'C']] = -4;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'Q']] = 2;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'G']] = -2;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'H']] = 0;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'I']] = -3;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'L']] = -3;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'K']] = 1;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'M']] = -2;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'F']] = -3;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'P']] = -1;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'S']] = 0;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'W']] = -3;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'Y']] = -2;
    b62[aa2int[(uint8_t) 'E']][aa2int[(uint8_t) 'V']] = -2;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'A']] = 0;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'R']] = -2;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'N']] = 0;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'D']] = -1;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'C']] = -3;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'Q']] = -2;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'E']] = -2;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'H']] = -2;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'I']] = -4;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'L']] = -4;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'K']] = -2;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'M']] = -3;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'F']] = -3;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'P']] = -2;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'S']] = 0;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'T']] = -2;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'W']] = -2;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'Y']] = -3;
    b62[aa2int[(uint8_t) 'G']][aa2int[(uint8_t) 'V']] = -3;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'A']] = -2;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'R']] = 0;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'N']] = 1;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'D']] = -1;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'C']] = -3;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'Q']] = 0;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'E']] = 0;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'G']] = -2;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'I']] = -3;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'L']] = -3;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'K']] = -1;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'M']] = -2;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'F']] = -1;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'P']] = -2;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'S']] = -1;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'T']] = -2;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'W']] = -2;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'Y']] = 2;
    b62[aa2int[(uint8_t) 'H']][aa2int[(uint8_t) 'V']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'A']] = -1;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'R']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'N']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'D']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'C']] = -1;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'Q']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'E']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'G']] = -4;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'H']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'L']] = 2;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'K']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'M']] = 1;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'F']] = 0;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'P']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'S']] = -2;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'W']] = -3;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'Y']] = -1;
    b62[aa2int[(uint8_t) 'I']][aa2int[(uint8_t) 'V']] = 3;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'A']] = -1;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'R']] = -2;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'N']] = -3;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'D']] = -4;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'C']] = -1;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'Q']] = -2;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'E']] = -3;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'G']] = -4;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'H']] = -3;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'I']] = 2;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'K']] = -2;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'M']] = 2;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'F']] = 0;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'P']] = -3;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'S']] = -2;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'W']] = -2;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'Y']] = -1;
    b62[aa2int[(uint8_t) 'L']][aa2int[(uint8_t) 'V']] = 1;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'A']] = -1;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'R']] = 2;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'N']] = 0;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'D']] = -1;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'C']] = -3;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'Q']] = 1;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'E']] = 1;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'G']] = -2;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'H']] = -1;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'I']] = -3;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'L']] = -2;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'M']] = -1;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'F']] = -3;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'P']] = -1;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'S']] = 0;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'W']] = -3;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'Y']] = -2;
    b62[aa2int[(uint8_t) 'K']][aa2int[(uint8_t) 'V']] = -2;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'A']] = -1;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'R']] = -1;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'N']] = -2;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'D']] = -3;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'C']] = -1;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'Q']] = 0;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'E']] = -2;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'G']] = -3;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'H']] = -2;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'I']] = 1;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'L']] = 2;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'K']] = -1;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'F']] = 0;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'P']] = -2;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'S']] = -1;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'W']] = -1;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'Y']] = -1;
    b62[aa2int[(uint8_t) 'M']][aa2int[(uint8_t) 'V']] = 1;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'A']] = -2;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'R']] = -3;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'N']] = -3;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'D']] = -3;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'C']] = -2;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'Q']] = -3;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'E']] = -3;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'G']] = -3;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'H']] = -1;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'I']] = 0;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'L']] = 0;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'K']] = -3;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'M']] = 0;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'P']] = -4;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'S']] = -2;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'T']] = -2;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'W']] = 1;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'Y']] = 3;
    b62[aa2int[(uint8_t) 'F']][aa2int[(uint8_t) 'V']] = -1;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'A']] = -1;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'R']] = -2;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'N']] = -2;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'D']] = -1;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'C']] = -3;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'Q']] = -1;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'E']] = -1;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'G']] = -2;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'H']] = -2;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'I']] = -3;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'L']] = -3;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'K']] = -1;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'M']] = -2;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'F']] = -4;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'S']] = -1;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'T']] = -1;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'W']] = -4;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'Y']] = -3;
    b62[aa2int[(uint8_t) 'P']][aa2int[(uint8_t) 'V']] = -2;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'A']] = 1;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'R']] = -1;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'N']] = 1;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'D']] = 0;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'C']] = -1;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'Q']] = 0;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'E']] = 0;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'G']] = 0;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'H']] = -1;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'I']] = -2;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'L']] = -2;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'K']] = 0;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'M']] = -1;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'F']] = -2;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'P']] = -1;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'T']] = 1;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'W']] = -3;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'Y']] = -2;
    b62[aa2int[(uint8_t) 'S']][aa2int[(uint8_t) 'V']] = -2;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'A']] = 0;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'R']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'N']] = 0;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'D']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'C']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'Q']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'E']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'G']] = -2;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'H']] = -2;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'I']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'L']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'K']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'M']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'F']] = -2;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'P']] = -1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'S']] = 1;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'W']] = -2;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'Y']] = -2;
    b62[aa2int[(uint8_t) 'T']][aa2int[(uint8_t) 'V']] = 0;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'A']] = -3;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'R']] = -3;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'N']] = -4;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'D']] = -4;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'C']] = -2;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'Q']] = -2;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'E']] = -3;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'G']] = -2;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'H']] = -2;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'I']] = -3;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'L']] = -2;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'K']] = -3;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'M']] = -1;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'F']] = 1;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'P']] = -4;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'S']] = -3;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'T']] = -2;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'Y']] = 2;
    b62[aa2int[(uint8_t) 'W']][aa2int[(uint8_t) 'V']] = -3;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'A']] = -2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'R']] = -2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'N']] = -2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'D']] = -3;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'C']] = -2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'Q']] = -1;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'E']] = -2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'G']] = -3;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'H']] = 2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'I']] = -1;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'L']] = -1;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'K']] = -2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'M']] = -1;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'F']] = 3;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'P']] = -3;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'S']] = -2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'T']] = -2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'W']] = 2;
    b62[aa2int[(uint8_t) 'Y']][aa2int[(uint8_t) 'V']] = -1;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'A']] = 0;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'R']] = -3;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'N']] = -3;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'D']] = -3;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'C']] = -1;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'Q']] = -2;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'E']] = -2;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'G']] = -3;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'H']] = -3;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'I']] = 3;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'L']] = 1;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'K']] = -2;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'M']] = 1;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'F']] = -1;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'P']] = -2;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'S']] = -2;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'T']] = 0;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'W']] = -3;
    b62[aa2int[(uint8_t) 'V']][aa2int[(uint8_t) 'Y']] = -1;

    // how about converting directly from the codon int representation to the AA int representation used by maxMAtches
    // (stop codons are just the max or min of that type..)
    // which in turn can be used to access blosum matrix directly
    // and could be synced with the aa int representation in bwt
    // then, to display sequences in debug mode, need to make a reverse map from int to char, but thats easy
    // those tables could all be in configkj and shared between threads, initialized in configkj constructor

    std::memset(codon2aa, '*', sizeof (codon2aa));

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("TCA")] = 'S';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable22) {
        codon2aa[codon_to_int("TCA")] = '*';
    } else {

    }

    codon2aa[codon_to_int("TCC")] = 'S';
    codon2aa[codon_to_int("TCG")] = 'S';
    codon2aa[codon_to_int("TCT")] = 'S';
    codon2aa[codon_to_int("TTC")] = 'F';
    codon2aa[codon_to_int("TTT")] = 'F';
    codon2aa[codon_to_int("TTA")] = 'L';
    codon2aa[codon_to_int("TTG")] = 'L';
    codon2aa[codon_to_int("TAC")] = 'Y';
    codon2aa[codon_to_int("TAT")] = 'Y';

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("TAA")] = '*';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable6 || mOptions->mTransSearchOptions->codonTable == codontable27) {
        codon2aa[codon_to_int("TAA")] = 'Q';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable29 ||
            mOptions->mTransSearchOptions->codonTable == codontable33 ||
            mOptions->mTransSearchOptions->codonTable == codontable14) {
        codon2aa[codon_to_int("TAA")] = 'Y';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable30) {
        codon2aa[codon_to_int("TAA")] = 'E';
    } else {

    }

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("TAG")] = '*';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable6 || mOptions->mTransSearchOptions->codonTable == codontable27) {
        codon2aa[codon_to_int("TAG")] = 'Q';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable16 ||
            mOptions->mTransSearchOptions->codonTable == codontable22) {
        codon2aa[codon_to_int("TAG")] = 'L';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable29) {
        codon2aa[codon_to_int("TAG")] = 'Y';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable30) {
        codon2aa[codon_to_int("TAG")] = 'E';
    } else {

    }

    codon2aa[codon_to_int("TGC")] = 'C';
    codon2aa[codon_to_int("TGT")] = 'C';

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("TGA")] = '*';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable10) {
        codon2aa[codon_to_int("TGA")] = 'C';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable2 ||
            mOptions->mTransSearchOptions->codonTable == codontable3 ||
            mOptions->mTransSearchOptions->codonTable == codontable4 ||
            mOptions->mTransSearchOptions->codonTable == codontable3 ||
            mOptions->mTransSearchOptions->codonTable == codontable5 ||
            mOptions->mTransSearchOptions->codonTable == codontable9 ||
            mOptions->mTransSearchOptions->codonTable == codontable14 ||
            mOptions->mTransSearchOptions->codonTable == codontable21 ||
            mOptions->mTransSearchOptions->codonTable == codontable24 ||
            mOptions->mTransSearchOptions->codonTable == codontable31 ||
            mOptions->mTransSearchOptions->codonTable == codontable33
            ) {
        codon2aa[codon_to_int("TGA")] = 'W';
    } else {
    }

    codon2aa[codon_to_int("TGG")] = 'W';

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("CTA")] = 'L';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable3) {
        codon2aa[codon_to_int("CTA")] = 'T';
    } else {

    }

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("CTC")] = 'L';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable3) {
        codon2aa[codon_to_int("CTC")] = 'T';
    } else {

    }


    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("CTG")] = 'L';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable3) {
        codon2aa[codon_to_int("CTG")] = 'T';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable12) {
        codon2aa[codon_to_int("CTG")] = 'S';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable26) {
        codon2aa[codon_to_int("CTG")] = 'A';
    } else {

    }

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("CTT")] = 'L';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable3) {
        codon2aa[codon_to_int("CTT")] = 'T';
    } else {

    }

    codon2aa[codon_to_int("CCA")] = 'P';
    codon2aa[codon_to_int("CAT")] = 'H';
    codon2aa[codon_to_int("CAA")] = 'Q';
    codon2aa[codon_to_int("CAG")] = 'Q';
    codon2aa[codon_to_int("CGA")] = 'R';
    codon2aa[codon_to_int("CGC")] = 'R';
    codon2aa[codon_to_int("CGG")] = 'R';
    codon2aa[codon_to_int("CGT")] = 'R';

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("ATA")] = 'I';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable2 ||
            mOptions->mTransSearchOptions->codonTable == codontable3 ||
            mOptions->mTransSearchOptions->codonTable == codontable5 ||
            mOptions->mTransSearchOptions->codonTable == codontable13 ||
            mOptions->mTransSearchOptions->codonTable == codontable21
            ) {
        codon2aa[codon_to_int("ATA")] = 'M';
    } else {

    }

    codon2aa[codon_to_int("ATC")] = 'I';
    codon2aa[codon_to_int("ATT")] = 'I';
    codon2aa[codon_to_int("ATG")] = 'M';
    codon2aa[codon_to_int("ACA")] = 'T';
    codon2aa[codon_to_int("ACC")] = 'T';
    codon2aa[codon_to_int("ACG")] = 'T';
    codon2aa[codon_to_int("ACT")] = 'T';
    codon2aa[codon_to_int("AAC")] = 'N';
    codon2aa[codon_to_int("AAT")] = 'N';

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("AAA")] = 'K';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable9 ||
            mOptions->mTransSearchOptions->codonTable == codontable14 ||
            mOptions->mTransSearchOptions->codonTable == codontable21) {
        codon2aa[codon_to_int("AAA")] = 'N';
    }


    codon2aa[codon_to_int("AAG")] = 'K';
    codon2aa[codon_to_int("AGC")] = 'S';
    codon2aa[codon_to_int("AGT")] = 'S';

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("AGA")] = 'R';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable2) {
        codon2aa[codon_to_int("AGA")] = '*';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable9 ||
            mOptions->mTransSearchOptions->codonTable == codontable5 ||
            mOptions->mTransSearchOptions->codonTable == codontable14 ||
            mOptions->mTransSearchOptions->codonTable == codontable21 ||
            mOptions->mTransSearchOptions->codonTable == codontable24 ||
            mOptions->mTransSearchOptions->codonTable == codontable33) {
        codon2aa[codon_to_int("AGA")] = 'S';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable13) {
        codon2aa[codon_to_int("AGA")] = 'G';
    } else {

    }

    if (mOptions->mTransSearchOptions->codonTable == codontable1) {
        codon2aa[codon_to_int("AGG")] = 'R';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable2) {
        codon2aa[codon_to_int("AGG")] = '*';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable5 ||
            mOptions->mTransSearchOptions->codonTable == codontable9 ||
            mOptions->mTransSearchOptions->codonTable == codontable14 ||
            mOptions->mTransSearchOptions->codonTable == codontable21 ||
            mOptions->mTransSearchOptions->codonTable == codontable24) {
        codon2aa[codon_to_int("AGG")] = 'S';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable13) {
        codon2aa[codon_to_int("AGG")] = 'G';
    } else if (mOptions->mTransSearchOptions->codonTable == codontable33) {
        codon2aa[codon_to_int("AGG")] = 'K';
    } else {

    }

    codon2aa[codon_to_int("CCC")] = 'P';
    codon2aa[codon_to_int("CCG")] = 'P';
    codon2aa[codon_to_int("CCT")] = 'P';
    codon2aa[codon_to_int("CAC")] = 'H';
    codon2aa[codon_to_int("GTA")] = 'V';
    codon2aa[codon_to_int("GTC")] = 'V';
    codon2aa[codon_to_int("GTG")] = 'V';
    codon2aa[codon_to_int("GTT")] = 'V';
    codon2aa[codon_to_int("GCA")] = 'A';
    codon2aa[codon_to_int("GCC")] = 'A';
    codon2aa[codon_to_int("GCG")] = 'A';
    codon2aa[codon_to_int("GCT")] = 'A';
    codon2aa[codon_to_int("GAC")] = 'D';
    codon2aa[codon_to_int("GAT")] = 'D';
    codon2aa[codon_to_int("GAA")] = 'E';
    codon2aa[codon_to_int("GAG")] = 'E';
    codon2aa[codon_to_int("GGA")] = 'G';
    codon2aa[codon_to_int("GGC")] = 'G';
    codon2aa[codon_to_int("GGG")] = 'G';
    codon2aa[codon_to_int("GGT")] = 'G';

    for (unsigned int i = 0; i <= 5; i++) {
        translations[i].reserve(5000);
    }
}

TransSearcher::~TransSearcher(){
}

Fragment *TransSearcher::getNextFragment(unsigned int min_score) {
    if (fragments.empty()) {
        return NULL;
    }
    auto it = fragments.begin();
    //if (mOptions->debug) std::cerr << "max fragment score/length = " << it->first << "\n";
    if (it->first < min_score) { //the highest scoring fragment in the sorted list is below threshold, then search stops
        return NULL;
    }
    Fragment *f = it->second;
    if (mOptions->debug) {
        std::stringstream ss;
        ss << "Fragment (protein) = " << f->seq << "; max fragment score/length = " << it->first;
        cCout(ss.str(), 'g');
    }
    fragments.erase(it);

    while (mOptions->mTransSearchOptions->comOptions.SEG && f != NULL && !f->SEGchecked) {
        std::string convertedseq = f->seq;
        for (size_t i = 0; i < convertedseq.length(); i++) {
            convertedseq[i] = AMINOACID_TO_NCBISTDAA[(int) convertedseq[i]];
        }
        BlastSeqLoc *seg_locs = NULL;
        SeqBufferSeg((Uint1 *) (convertedseq.data()), (Int4) convertedseq.length(), 0, tbwtfmiDB->blast_seg_params, &seg_locs);
        if (seg_locs) { // SEG found region(s)
            BlastSeqLoc *curr_loc = seg_locs;
            size_t start = 0; //start of non-SEGged piece
            do {
                size_t length = curr_loc->ssr->left - start;
                if (mOptions->debug) {
                    std::stringstream ss;
                    ss << "SEG region: " << curr_loc->ssr->left << " - " << curr_loc->ssr->right << " = " << f->seq.substr(curr_loc->ssr->left, curr_loc->ssr->right - curr_loc->ssr->left + 1);
                    cCout(ss.str(), 'g');
                }
                //if (length > mOptions->mTransSearchOptions->comOptions.minFragLength) {
                if (length > min_match_length) {
                    if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                        unsigned int score = calcScore(f->seq, start, length, 0);
                        if (score >= mOptions->mTransSearchOptions->comOptions.minScore) {
                            fragments.emplace(score, new Fragment(f->seq.substr(start, length), true));
                        }
                    } else {
                        fragments.emplace(length, new Fragment(f->seq.substr(start, length), true));
                    }
                }
                start = curr_loc->ssr->right + 1;
            } while ((curr_loc = curr_loc->next) != NULL);
            size_t len_last_piece = f->seq.length() - start;
            //if (len_last_piece > mOptions->mTransSearchOptions->comOptions.minFragLength) {
            if (len_last_piece > min_match_length) {
                if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                    unsigned int score = calcScore(f->seq, start, len_last_piece, 0);
                    if (score >= mOptions->mTransSearchOptions->comOptions.minScore) {
                        fragments.emplace(score, new Fragment(f->seq.substr(start, len_last_piece), true));
                    }
                } else {
                    fragments.emplace(len_last_piece, new Fragment(f->seq.substr(start, len_last_piece), true));
                }
            }

            BlastSeqLocFree(seg_locs);
            delete f;
            f = NULL;
            if (!fragments.empty()) {
                it = fragments.begin();
                if (it->first >= min_score) {
                    f = it->second;
                    fragments.erase(it);
                    // next iteration of while loop
                }
            }
        } else { // no SEG regions found
            return f;
        }
    }

    return f;
}

void TransSearcher::getAllFragmentsBits(const std::string& line) {

    for (unsigned int i = 0; i <= 2; i++) {
        translations[i].clear();
    }
    const char *c = line.c_str();
    for (size_t count = 0; count < line.length() - 2; count++) {
        char aa = codon2aa[codon_to_int(c++)];
        if (aa == '*') {
            size_t index = count % 3;
            // finished one of the translations, so add it to fragments
            //if (translations[index].length() >= mOptions->mTransSearchOptions->comOptions.minFragLength) {
            if (translations[index].length() >= min_match_length) {
                if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                    unsigned int score = calcScore(translations[index]);
                    if (score >= mOptions->mTransSearchOptions->comOptions.minScore)
                        fragments.emplace(score, new Fragment(translations[index]));
                } else {
                    fragments.emplace(translations[index].length(), new Fragment(translations[index]));
                }
            }
            translations[index].clear();
        } else {
            translations[count % 3].push_back(aa);
        }
    }

    for (unsigned int i = 0; i <= 2; i++) {
        //add remaining stuff to fragments
        //if (translations[i].length() >= mOptions->mTransSearchOptions->comOptions.minFragLength) {
        if (translations[i].length() >= min_match_length) {
            if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                unsigned int score = calcScore(translations[i]);
                if (score >= mOptions->mTransSearchOptions->comOptions.minScore)
                    fragments.emplace(score, new Fragment(translations[i]));
            } else {
                fragments.emplace(translations[i].length(), new Fragment(translations[i]));
            }
        }
        translations[i].clear();
    }

    //cerr << "\n";
    for (int count = (int) line.length() - 2; count >= 0; count--) { // count needs to be an int here and not size_t
        //	printf("=%i.",revcomp_codon_to_int(codon));
        char aa = codon2aa[revcomp_codon_to_int(c--)];
        if (aa == '*') {
            size_t index = count % 3;
            // finished one of the translations, so add it to fragments
            //if (translations[index].length() >= mOptions->mTransSearchOptions->comOptions.minFragLength) {
            if (translations[index].length() >= min_match_length) {
                if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                    unsigned int score = calcScore(translations[index]);
                    if (score >= mOptions->mTransSearchOptions->comOptions.minScore)
                        fragments.emplace(score, new Fragment(translations[index]));
                } else {
                    fragments.emplace(translations[index].length(), new Fragment(translations[index]));
                }
            }
            translations[index].clear();
        } else {
            translations[count % 3].push_back(aa);
        }
    }

    for (unsigned int i = 0; i <= 2; i++) {
        //add remaining stuff to fragments
        //if (translations[i].length() >= mOptions->mTransSearchOptions->comOptions.minFragLength) {
        if (translations[i].length() >= min_match_length) {
            if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                unsigned int score = calcScore(translations[i]);
                if (score >= mOptions->mTransSearchOptions->comOptions.minScore)
                    fragments.emplace(score, new Fragment(translations[i]));
            } else {
                fragments.emplace(translations[i].length(), new Fragment(translations[i]));
            }
        }
    }
}

void TransSearcher::getLongestFragmentsBits(const std::string& line) {

    unsigned int min_len_cutoff = 0;

    if (line.length() % 3 == 2) {
        min_len_cutoff = floor(line.length() / 3);
    } else {
        min_len_cutoff = floor(line.length() / 3) - 1;
    }

    min_len_cutoff = max(min_len_cutoff, min_match_length);
    //min_len_cutoff = max(min_len_cutoff, mOptions->mTransSearchOptions->comOptions.minFragLength);

    for (unsigned int i = 0; i <= 2; i++) {
        translations[i].clear();
    }
    const char *c = line.c_str();
    for (size_t count = 0; count < line.length() - 2; count++) {
        char aa = codon2aa[codon_to_int(c++)];
        if (aa == '*') {
            size_t index = count % 3;
            // finished one of the translations, so add it to fragments
            if (translations[index].length() >= min_len_cutoff) {
                if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                    unsigned int score = calcScore(translations[index]);

                    if (score >= mOptions->mTransSearchOptions->comOptions.minScore)
                        fragments.emplace(score, new Fragment(translations[index]));
                } else {
                    fragments.emplace(translations[index].length(), new Fragment(translations[index]));
                }
            }
            translations[index].clear();
        } else {
            translations[count % 3].push_back(aa);
        }
    }
    for (unsigned int i = 0; i <= 2; i++) {
        //add remaining stuff to fragments
        if (translations[i].length() >= min_len_cutoff) {
            if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                unsigned int score = calcScore(translations[i]);
                if (score >= mOptions->mTransSearchOptions->comOptions.minScore)
                    fragments.emplace(score, new Fragment(translations[i]));
            } else {
                fragments.emplace(translations[i].length(), new Fragment(translations[i]));
            }
        }
        translations[i].clear();
    }
    //cerr << "\n";
    for (int count = (int) line.length() - 2; count >= 0; count--) { // count needs to be an int here and not size_t
        //	printf("=%i.",revcomp_codon_to_int(codon));
        char aa = codon2aa[revcomp_codon_to_int(c--)];
        if (aa == '*') {
            size_t index = count % 3;
            // finished one of the translations, so add it to fragments
            if (translations[index].length() >= min_len_cutoff) {
                if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                    unsigned int score = calcScore(translations[index]);
                    if (score >= mOptions->mTransSearchOptions->comOptions.minScore)
                        fragments.emplace(score, new Fragment(translations[index]));
                } else {
                    fragments.emplace(translations[index].length(), new Fragment(translations[index]));
                }
            }
            translations[index].clear();
        } else {
            translations[count % 3].push_back(aa);
        }
    }
    for (unsigned int i = 0; i <= 2; i++) {
        //add remaining stuff to fragments
        if (translations[i].length() >= min_len_cutoff) {
            if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
                unsigned int score = calcScore(translations[i]);
                if (score >= mOptions->mTransSearchOptions->comOptions.minScore)
                    fragments.emplace(score, new Fragment(translations[i]));
            } else {
                fragments.emplace(translations[i].length(), new Fragment(translations[i]));
            }
        }
    }
}

void TransSearcher::addAllMismatchVariantsAtPosSI(const Fragment *f, unsigned int pos, size_t erase_pos = std::string::npos, SI *si = NULL) {

    assert(mOptions->mTransSearchOptions->comOptions.mode == GREEDY);
    assert(pos < erase_pos);
    assert(f->num_mm == 0 || pos < f->pos_lastmm);

    std::string fragment = f->seq; // make a copy to modify the sequence at pos
    assert(fragment.length() >= min_match_length);
    //assert(fragment.length() >= mOptions->mTransSearchOptions->comOptions.minFragLength);
    char origchar = fragment[pos];
    assert(blosum_subst.count(origchar) > 0);

    if (erase_pos != std::string::npos && erase_pos < fragment.length()) {
        if (mOptions->debug) {
            std::stringstream ss;
            ss << "Deleting from position " << erase_pos;
            cCout(ss.str(), 'g');
        }
        fragment.erase(erase_pos);
    }

    //calc score for whole sequence, so we can substract the diff for each substitution
    unsigned int score = calcScore(fragment, f->diff) - blosum62diag[aa2int[(uint8_t) origchar]];
    IndexType siarray[2], siarrayupd[2];
    siarray[0] = si->start;
    siarray[1] = si->start + (IndexType) si->len;

    for (auto itv : blosum_subst.at(origchar)) {
        // we know the difference between score of original aa and substitution score, this
        // has to be subtracted when summing over all positions later
        // so we add this difference to the fragment
        int score_after_subst = score + b62[aa2int[(uint8_t) origchar]][aa2int[(uint8_t) itv]];
        if (score_after_subst >= (int) best_match_score && score_after_subst >= (int) mOptions->mTransSearchOptions->comOptions.minScore) {
            if (UpdateSI(tbwtfmiDB->fmi, tbwtfmiDB->astruct->trans[(size_t) itv], siarray, siarrayupd, 'p') != 0) {
                fragment[pos] = itv;
                int diff = b62[aa2int[(uint8_t) origchar]][aa2int[(uint8_t) itv]] - blosum62diag[aa2int[(uint8_t) itv]];
                if (mOptions->debug) {
                    std::stringstream ss;
                    ss << "Adding fragment   " << fragment << " with mismatch at pos " << pos << " ,diff " << f->diff + diff << ", max score " << score_after_subst;
                    cCout(ss.str(), 'g');
                }
                fragments.emplace(score_after_subst, new Fragment(fragment, f->num_mm + 1, pos, f->diff + diff, siarrayupd[0], siarrayupd[1], si->ql + 1));
            } else if (mOptions->debug) {
                fragment[pos] = itv;
                std::stringstream ss;
                ss << "Skipping fragment " << fragment << " mismatch at pos " << pos << ", because " << itv << " is not a valid extension";
                cCout(ss.str(), 'g');
            }
        } else {
            if (mOptions->debug) {
                fragment[pos] = itv;
                std::stringstream ss;
                ss << "Skipping fragment " << fragment << " and following fragments, because score is too low: " << score_after_subst << " < " << std::max(best_match_score, mOptions->mTransSearchOptions->comOptions.minScore);
                cCout(ss.str(), 'g');
            }
            break;
        }
    }
}

unsigned int TransSearcher::calcScore(const std::string &s, size_t start, size_t len, int diff) {
    int score = 0;
    for (size_t i = start; i < start + len; ++i) {
        score += blosum62diag[aa2int[(uint8_t) s[i]]];
    }
    score += diff;
    return score > 0 ? score : 0;
}

unsigned int TransSearcher::calcScore(const std::string &s, int diff) {
    int score = 0;
    for (size_t i = 0; i < s.length(); ++i) {
        score += blosum62diag[aa2int[(uint8_t) s[i]]];
    }
    score += diff;
    return score > 0 ? score : 0;
}

unsigned int TransSearcher::calcScore(const std::string &s) {
    unsigned int score = 0;
    for (size_t i = 0; i < s.length(); ++i) {
        score += blosum62diag[aa2int[(uint8_t) s[i]]];
    }
    return score;
}

void TransSearcher::eval_match_scores(SI *si, Fragment *frag) {

    if (!si)
        return;

    // eval the remaining same-length and shorter matches
    if (si->samelen)
        eval_match_scores(si->samelen, frag);
    //if (si->next && si->next->ql >= (int) mOptions->mTransSearchOptions->comOptions.minFragLength)
    if (si->next && si->next->ql >= (int) min_match_length)
        eval_match_scores(si->next, frag);
    else if (si->next)
        recursive_free_SI(si->next);

    unsigned int score = calcScore(frag->seq, si->qi, si->ql, frag->diff);

    if (mOptions->debug) {
        std::stringstream ss;
        ss << "Match " << frag->seq.substr(si->qi, si->ql) << " (length=" << (unsigned int) si->ql << " score=" << score << " num_mm=" << frag->num_mm << ")";
        cCout(ss.str(), 'g');
    }

    if (score < mOptions->mTransSearchOptions->comOptions.minScore) {
        free(si);
        si = NULL;
        return;
    }

    if (score > best_match_score) {
        for (auto itm : best_matches_SI) {
            //recursive_free_SI(itm);
            free(itm);
        }
        best_matches_SI.clear();
        best_matches_SI.push_back(si);
        best_match_score = score;
        // match_str.append(frag->seq.substr(si->qi, si->ql));
        // match_str.append(";");

        if (mOptions->verbose) {
            best_matches.clear();
            best_matches.push_back(frag->seq.substr(si->qi, si->ql));
        }
    } else if (score == best_match_score && best_matches_SI.size() < mOptions->mTransSearchOptions->comOptions.max_matches_SI) {
        best_matches_SI.push_back(si);
        // match_str.append(frag->seq.substr(si->qi, si->ql));
        // match_str.append(";");
        if (mOptions->verbose)
            best_matches.push_back(frag->seq.substr(si->qi, si->ql));
    } else {
        free(si);
        si = NULL;
    }
}

void TransSearcher::clearFragments() {
    for(auto it = fragments.begin(); it != fragments.end(); ++it){
        delete it->second;
    }
    fragments.clear();
}

void TransSearcher::clearMatchedIds(){
    for(auto it = match_ids.begin(); it != match_ids.end(); ++it) {
        if(*it){
            delete *it;
        }
    }
    match_ids.clear();
}

inline uint8_t TransSearcher::codon_to_int(const char *codon) {
    return (uint8_t) (nuc2int[(uint8_t) codon[0]] << 4 | nuc2int[(uint8_t) codon[1]] << 2 | nuc2int[(uint8_t) codon[2]]);
}

inline uint8_t TransSearcher::revcomp_codon_to_int(const char *codon) {
    return (uint8_t) (compnuc2int[(uint8_t) codon[2]] << 4 | compnuc2int[(uint8_t) codon[1]] << 2 | compnuc2int[(uint8_t) codon[0]]);
}

void TransSearcher::classify_greedyblosum() {
    best_matches_SI.clear();
    best_matches.clear();
    best_match_score = 0;

    while (1) {
        Fragment *t = getNextFragment(best_match_score);
        if (!t) break;
        const std::string fragment = t->seq;
        const size_t length = fragment.length();
        const unsigned int num_mm = t->num_mm;

        if (mOptions->debug) {
            std::stringstream ss;
            ss << "Searching fragment " << fragment << " (" << length << "," << num_mm << "," << t->diff << ")";
            cCout(ss.str(), 'g');
        }
        char *seq = new char[length + 1];
        std::strcpy(seq, fragment.c_str());
        //stopped here.
        translate2numbers((uchar *) seq, (unsigned int) length, tbwtfmiDB->astruct);
        SI *si = NULL;
        if (num_mm > 0) {
            if (num_mm == mOptions->mTransSearchOptions->comOptions.misMatches) { //after last mm has been done, we need to have at least reached the min_length
                //si = maxMatches_withStart(tbwtfmiDB->fmi, seq, (unsigned int) length, mOptions->mTransSearchOptions->comOptions.minFragLength, 1, t->si0, t->si1, t->matchlen, 'p');
                si = maxMatches_withStart(tbwtfmiDB->fmi, seq, (unsigned int) length, min_match_length, 1, t->si0, t->si1, t->matchlen, 'p');
            } else {
                si = maxMatches_withStart(tbwtfmiDB->fmi, seq, (unsigned int) length, t->matchlen, 1, t->si0, t->si1, t->matchlen, 'p');
            }
        } else {
            si = maxMatches(tbwtfmiDB->fmi, seq, (unsigned int) length, mOptions->mTransSearchOptions->comOptions.seedLength, 0, 'p'); //initial matches
        }
        if (!si) { // no match for this fragment
            if (mOptions->debug) {
                std::stringstream ss;
                ss << "No match for this fragment.";
                cCout(ss.str(), 'g');
            }
            delete[] seq;
            delete t;
            continue; // continue with the next fragment
        }
        if (mOptions->debug) {
            std::stringstream ss;
            ss << "Longest match has length " << (unsigned int) si->ql;
            cCout(ss.str(), 'g');
        }

        if (mOptions->mTransSearchOptions->comOptions.misMatches > 0 && num_mm < mOptions->mTransSearchOptions->comOptions.misMatches) {
            SI *si_it = si;
            while (si_it) {
                unsigned int match_right_end = si_it->qi + si_it->ql - 1;
                if (num_mm > 0)
                    assert(match_right_end == length - 1); // greedy matches end always at the end
                if (mOptions->debug) {
                    std::stringstream ss;
                    ss << "Match from " << si_it->qi << " to " << match_right_end << ": " << fragment.substr(si_it->qi, match_right_end - si_it->qi + 1) << " (" << si_it->ql << ")";
                    cCout(ss.str(), 'g');
                }
                //if (si_it->qi > 0 && match_right_end + 1 >= mOptions->mTransSearchOptions->comOptions.minFragLength) {
                if (si_it->qi > 0 && match_right_end + 1 >= min_match_length) {
                    //1. match must end before beginning of fragment, i.e. it is extendable
                    //2. remaining fragment, from zero to end of current match, must be longer than minimum length of accepted matches
                    const size_t erase_pos = (match_right_end < length - 1) ? match_right_end + 1 : std::string::npos;
                    addAllMismatchVariantsAtPosSI(t, (unsigned int) (si_it->qi - 1), erase_pos, si_it);
                }
                si_it = si_it->samelen ? si_it->samelen : si_it->next;
            }
        }

        //if ((unsigned int) si->ql < mOptions->mTransSearchOptions->comOptions.minFragLength) { // match was too short
        if ((unsigned int) si->ql < min_match_length) { // match was too short
            if (mOptions->debug) {
                std::stringstream ss;
                ss << "Match of length " << si->ql << " is too short";
                cCout(ss.str(), 'g');
            }
            delete[] seq;
            delete t;
            recursive_free_SI(si);
            continue; // continue with the next fragment
        }

        eval_match_scores(si, t);

        delete[] seq;
        delete t;

    } // end current fragment

    if (best_matches_SI.empty()) {
        return;
    }

    if (mOptions->mTransSearchOptions->comOptions.useEvalue) {
        //calc e-value and only return match if > cutoff

        double bitscore = (LAMBDA * best_match_score - LN_K) / LN_2;
        double Evalue = tbwtfmiDB->db_length * query_len * pow(2, -1 * bitscore);
        if (mOptions->debug) {
            std::stringstream ss;
            ss << "E-value = " << Evalue;
            cCout(ss.str(), 'g');
        }

        if (Evalue > mOptions->mTransSearchOptions->comOptions.minEvalue) {
            for (auto itm : best_matches_SI) {
                free(itm);
                //recursive_free_SI(itm);
            }
            return;
        }
    }

    match_ids.clear();

    for (auto itm : best_matches_SI) {
        ids_from_SI(itm);
    }
    for (auto itm : best_matches_SI) {
        //recursive_free_SI(itm); // why not use this one;
         free(itm);
    }
}

void TransSearcher::classify_length() {
    unsigned int longest_match_length = 0;
    longest_matches_SI.clear();
    longest_fragments.clear();

    while (1) {
        Fragment *t = getNextFragment(longest_match_length);
        if (!t) break; // searched all fragments that are longer than best match length
        const std::string fragment = t->seq;
        const unsigned int length = (unsigned int) fragment.length();

        if (mOptions->debug) {
            std::stringstream ss;
            ss << "Searching fragment (protein) " << fragment << " (" << length << ")";
            cCout(ss.str(), 'g');
        }
        char *seq = new char[length + 1];
        std::strcpy(seq, fragment.c_str());

        translate2numbers((uchar *) seq, length, tbwtfmiDB->astruct);
        //use longest_match_length here too:
        //SI *si = maxMatches(tbwtfmiDB->fmi, seq, length, std::max(mOptions->mTransSearchOptions->comOptions.minFragLength, longest_match_length), 1, 'p');
        SI *si = maxMatches(tbwtfmiDB->fmi, seq, length, std::max(min_match_length, longest_match_length), 1, 'p');

        if (!si) { // no match for this fragment
            if (mOptions->debug) {
                std::stringstream ss;
                ss << "No match for this fragment.";
                cCout(ss.str(), 'g');
            }
            delete[] seq;
            delete t;
            continue; // continue with the next fragment
        }

        // just get length here and save si when it is longest
        if (mOptions->debug) {
            std::stringstream ss;
            ss << "Longest match is length " << (unsigned int) si->ql;
            cCout(ss.str(), 'g');
        }
        if ((unsigned int) si->ql > longest_match_length) {
            for (auto itm : longest_matches_SI)
                recursive_free_SI(itm);
            longest_matches_SI.clear();
            longest_matches_SI.push_back(si);
            longest_match_length = (unsigned int) si->ql;
            // match_str.append(fragment.substr(si->qi, si->ql));
            // match_str.append(";");
            if (mOptions->verbose) {
                longest_fragments.clear();
                longest_fragments.push_back(fragment.substr(si->qi, si->ql));
            }
        } else if ((unsigned int) si->ql == longest_match_length) {
            longest_matches_SI.push_back(si);
            // match_str.append(fragment.substr(si->qi, si->ql));
            // match_str.append(";");
            if (mOptions->verbose)
                longest_fragments.push_back(fragment.substr(si->qi, si->ql));
        } else {
            recursive_free_SI(si);
            si = NULL;
        }
        delete[] seq;
        delete t;

    } // end current fragment

    if (longest_matches_SI.empty()) {
        return;
    }
    match_ids.clear();
    for (auto itm : longest_matches_SI) {
        ids_from_SI_recursive(itm);
    }
    for (auto itm : longest_matches_SI) {
        recursive_free_SI(itm);
    }
}

void TransSearcher::postProcess() {
}

std::set<char *>& TransSearcher::transSearchWuKong(Read* & item) {
    //match_str.clear();
    min_match_length = static_cast<uint>(mOptions->mTransSearchOptions->comOptions.lenper * item->length() / 3);
    clearFragments();
    //clearMatchedIds();
    match_ids.clear();
    query_len = static_cast<double> (item->length()) / 3.0;
    if (mOptions->debug) {
        std::stringstream ss;
        ss << "Getting fragments for read: " << item->mName << "\t" << item->mSeq.mStr;
        cCout(ss.str(), 'g');
    }
    getAllFragmentsBits(item->mSeq.mStr);

    if (mOptions->debug) {
        std::stringstream ss;
        ss << fragments.size() << " fragments found in the read";
        cCout(ss.str(), 'g');
    }

    if (mOptions->mTransSearchOptions->comOptions.mode == MEM) {
        classify_length();
    } else if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
        classify_greedyblosum();
    } else { // this should not happen
        assert(false);
    }

    clearFragments();
    return match_ids;
}

std::set<char *>& TransSearcher::transSearchWuKong(Read* & item1, Read* & item2) {
    //match_str.clear();
    min_match_length = static_cast<uint>(mOptions->mTransSearchOptions->comOptions.lenper * std::min(item1->length(), item2->length()) / 3);
    clearFragments();
    //clearMatchedIds();
    match_ids.clear();
    query_len = static_cast<double> (item1->length()) / 3.0;
    if (mOptions->debug) {
        std::stringstream ss;
        ss << "Getting fragments for read1: " << item1->mName << "\t" << item1->mSeq.mStr;
        cCout(ss.str(), 'g');
    }
    getAllFragmentsBits(item1->mSeq.mStr);

    query_len += static_cast<double> (item2->length()) / 3.0;
    if (mOptions->debug) {
        std::stringstream ss;
        ss << "Getting fragments for read2: " << item2->mName << "\t" << item2->mSeq.mStr;
        cCout(ss.str(), 'g');
    }
    getAllFragmentsBits(item2->mSeq.mStr);

    if (mOptions->debug) {
        std::stringstream ss;
        ss << fragments.size() << " fragments found in the read.";
        cCout(ss.str(), 'g');
    }

    if (mOptions->mTransSearchOptions->comOptions.mode == MEM) {
        classify_length();
    } else if (mOptions->mTransSearchOptions->comOptions.mode == GREEDY) {
        classify_greedyblosum();
    } else { // this should not happen
        assert(false);
    }
    clearFragments();
    return match_ids;
}

void TransSearcher::ids_from_SI(SI *si) {
    IndexType k, pos;
    int iseq;
    for (k = si->start; k < si->start + si->len; ++k) {
        if (match_ids.size() > mOptions->mTransSearchOptions->comOptions.max_match_ids) {
            break;
        }
        get_suffix(tbwtfmiDB->fmi, tbwtfmiDB->bwt->s, k, &iseq, &pos, 'p');
        match_ids.insert(tbwtfmiDB->bwt->s->ids[iseq]);
    }
}

void TransSearcher::ids_from_SI_recursive(SI *si) {
    SI *si_it = si;
    while (si_it) {
        IndexType k, pos;
        int iseq;
        for (k = si_it->start; k < si_it->start + si_it->len; ++k) {
            if (match_ids.size() > mOptions->mTransSearchOptions->comOptions.max_match_ids) {
                break;
            }
            get_suffix(tbwtfmiDB->fmi, tbwtfmiDB->bwt->s, k, &iseq, &pos, 'p');
            match_ids.insert(tbwtfmiDB->bwt->s->ids[iseq]);
        } // end for
        si_it = si_it->samelen;
    } // end while all SI with same length
}