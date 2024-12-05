#include "dnasearcher.h"

DNASearcher::DNASearcher(Options * & opt, BwtFmiDB* & bwtfmiDB) {
    mOptions = opt;
    mBwtfmiDB = bwtfmiDB;
    std::memset(nuc2int, std::numeric_limits<uint8_t>::max(), sizeof(nuc2int));
    nuc2int['A'] = nuc2int['a'] = 0;
    nuc2int['C'] = nuc2int['c'] = 1;
    nuc2int['T'] = nuc2int['g'] = 2;
    nuc2int['G'] = nuc2int['t'] = 3;
    nuc2int['U'] = nuc2int['u'] = 3;

    std::memset(compnuc2int, std::numeric_limits<uint8_t>::max(), sizeof (compnuc2int));
    compnuc2int['A'] = compnuc2int['a'] = 3;
    compnuc2int['C'] = compnuc2int['c'] = 2;
    compnuc2int['G'] = compnuc2int['g'] = 1;
    compnuc2int['T'] = compnuc2int['t'] = 0;
    compnuc2int['U'] = compnuc2int['u'] = 0;

    blosum_subst = {
        {'A', {'C', 'G', 'T'}},
        {'C', {'A', 'G', 'T'}},
        {'G', {'A', 'C', 'T'}},
        {'T', {'A', 'C', 'G'}}
    };

    blosum62diag[nuc2int['A']] = 1;
    blosum62diag[nuc2int['C']] = 1;
    blosum62diag[nuc2int['G']] = 1;
    blosum62diag[nuc2int['T']] = 1;

    b62[nuc2int[(uint8_t) 'A']][nuc2int[(uint8_t) 'C']] = 0;
    b62[nuc2int[(uint8_t) 'A']][nuc2int[(uint8_t) 'T']] = 0;
    b62[nuc2int[(uint8_t) 'A']][nuc2int[(uint8_t) 'G']] = 1;

    b62[nuc2int[(uint8_t) 'C']][nuc2int[(uint8_t) 'A']] = 0;
    b62[nuc2int[(uint8_t) 'C']][nuc2int[(uint8_t) 'G']] = 0;
    b62[nuc2int[(uint8_t) 'C']][nuc2int[(uint8_t) 'T']] = 1;

    b62[nuc2int[(uint8_t) 'G']][nuc2int[(uint8_t) 'A']] = 1;
    b62[nuc2int[(uint8_t) 'G']][nuc2int[(uint8_t) 'C']] = 0;
    b62[nuc2int[(uint8_t) 'G']][nuc2int[(uint8_t) 'T']] = 0;

    b62[nuc2int[(uint8_t) 'T']][nuc2int[(uint8_t) 'A']] = 0;
    b62[nuc2int[(uint8_t) 'T']][nuc2int[(uint8_t) 'C']] = 1;
    b62[nuc2int[(uint8_t) 'T']][nuc2int[(uint8_t) 'G']] = 0;
}

DNASearcher::~DNASearcher(){
    
}

void DNASearcher::clearFragments() {
    for(auto it = fragments.begin(); it != fragments.end(); ++it){
        delete it->second;
    }
    fragments.clear();
}

void DNASearcher::clearMatchedIds(){
    for(auto it = match_ids.begin(); it != match_ids.end(); ++it) {
        if(*it){
            delete *it;
        }
    }
    match_ids.clear();
}

void DNASearcher::getAllFragments(Read * & item) {
    Read* rr1 = nullptr;

    // rr1 = item->reverseComplement();
    // fragments.emplace(item->length(), new Fragment(item->mSeq.mStr));
    // fragments.emplace(rr1->length(), new Fragment(rr1->mSeq.mStr));
    if (mOptions->mDNASearchOptions->mode == dGREEDY) {
        uint32 score = calcScore(item->mSeq.mStr);
        if (score >= mOptions->mDNASearchOptions->minScore) {
            fragments.emplace(score, new Fragment(item->mSeq.mStr));
        }
        rr1 = item->reverseComplement();
        score = calcScore(rr1->mSeq.mStr);
        if (score >= mOptions->mDNASearchOptions->minScore) {
            fragments.emplace(score, new Fragment(rr1->mSeq.mStr));
        }
    } else {
        rr1 = item->reverseComplement();
        fragments.emplace(item->length(), new Fragment(item->mSeq.mStr));
        fragments.emplace(rr1->length(), new Fragment(rr1->mSeq.mStr));
    }
    if(rr1) delete rr1; rr1 = nullptr;
}

std::set<char *>& DNASearcher::dnaSearchWuNeng(Read* & item) {
    clearFragments();
    //clearMatchedIds();
    match_ids.clear();
    query_len = static_cast<double>(item->length());
    getAllFragments(item);
    if (mOptions->mDNASearchOptions->mode == dMEM) {
        classify_length();
    } else if (mOptions->mDNASearchOptions->mode == dGREEDY) {
        classify_greedyblosum();
    } else { // this should not happen
        assert(false);
    }
    clearFragments();
    if (!match_ids.empty()) {
        //postProcess();
    }
    return (match_ids);
}

std::set<char *>& DNASearcher::dnaSearchWuNeng(Read* & item1, Read* & item2) {
    clearFragments();
    //clearMatchedIds();
    match_ids.clear();
    readName = item1->mName;
    query_len = static_cast<double>(item1->length());
    getAllFragments(item1);
    query_len += static_cast<double> (item2->length());
    getAllFragments(item2);

    if (mOptions->mDNASearchOptions->mode == dMEM) {
        classify_length();
        //cCout("length", 'r');
    } else if (mOptions->mDNASearchOptions->mode == dGREEDY) {
        classify_greedyblosum();
    } else { // this should not happen
        assert(false);
    }

    clearFragments();

    if (!match_ids.empty()) {
        //postProcess();
    }
    return (match_ids);
}

void DNASearcher::postProcess() {
    if (mOptions->verbose) {
        ss.str("");
        ss << readName << "\t";

        for (const auto it : match_ids) {
            if (it == *(match_ids.rbegin())) {
                ss << it;
            } else {
                ss << it << ";";
            }
        }
        cCout(ss.str(), 'g');
    }
}

void DNASearcher::classify_length() {
    uint longest_match_length = 0;
    longest_matches_SI.clear();
    longest_fragments.clear();
    //if (mOptions->debug) std::cerr << "debug: classify_length: a" << "\n";
    while (1) {
        Fragment* t = getNextFragment(longest_match_length);//need to delete t;
        if (!t) break;
        const std::string fragment = t->seq;
        const uint length = (uint) fragment.length();

        if (mOptions->debug) {
            std::stringstream ss;
            ss << "Searching DNA fragment " << fragment << " (" << length << ")" << "\n";
            cCout(ss.str(), 'y');
        }

        char *seq = new char[length + 1];//need to delete seq
        std::strcpy(seq, fragment.c_str());
        translate2numbers((uchar*) seq, length, mBwtfmiDB->astruct);

        //need to delete si;
        SI* si = maxMatches(mBwtfmiDB->fmi, seq, length, std::max(mOptions->mDNASearchOptions->minFragLength, longest_match_length), 1);

        if (!si) {
            if (mOptions->debug) {
                std::stringstream ss;
                ss << "No match for this fragment.";
                cCout(ss.str(), 'y');
            }
            delete[] seq;
            delete t;
            continue;
        } 

        //if (mOptions->debug) std::cerr << "Longest match (DNA) is length " << (uint) si->ql << "\n";

        if ((uint) si->ql > longest_match_length) {
            for (auto itm : longest_matches_SI) recursive_free_SI(itm);
            longest_matches_SI.clear();
            longest_matches_SI.push_back(si);
            longest_match_length = (uint) si->ql;
            if (mOptions->verbose) {
                longest_fragments.clear();
                longest_fragments.push_back(fragment.substr(si->qi, si->ql));
            }
        } else if ((uint) si->ql == longest_match_length) {
            longest_matches_SI.push_back(si);
            if (mOptions->verbose) longest_fragments.push_back(fragment.substr(si->qi, si->ql));
        } else {
            recursive_free_SI(si);
            si = NULL;
        }

        // if (mOptions->debug) std::cerr << "Longest match (DNA) finished " << (uint) si->ql << "\n";

        delete[] seq;
        delete t;

        // if (mOptions->debug) std::cerr << "Longest match (DNA) finished3 " << (uint) si->ql << "\n";
    }

    //if (mOptions->debug) std::cerr << "debug: classify_length: b" << "\n";
    if (longest_matches_SI.empty()) {
        return;
    }

    //if (mOptions->debug) std::cerr << "debug: classify_length: c" << "\n";
    match_ids.clear();

    //if (mOptions->debug) std::cerr << "debug: classify_length: d" << "\n";
    for (auto itm : longest_matches_SI) {
        ids_from_SI_recursive(itm);
    }

    // if (mOptions->debug) std::cerr << "debug: classify_length: e" << "\n";
    for (auto itm : longest_matches_SI) {
        recursive_free_SI(itm);
    }
    // if (mOptions->debug) std::cerr << "debug: classify_length: f" << "\n";
}

void DNASearcher::classify_greedyblosum() {
    best_matches_SI.clear();
    best_matches.clear();
    best_match_score = 0;

    while (1) {
        Fragment *t = getNextFragment(best_match_score);
        if (!t) break;

        const std::string fragment = t->seq;
        const size_t length = fragment.length();
        const uint32 num_mm = t->num_mm;

        if (mOptions->debug) {
            std::stringstream ss;
            ss << "Searching fragment (DNA) " << fragment << " (" << length << "," << num_mm << "," << t->diff << ")";
            cCout(ss.str(), 'y');
        }

        char *seq = new char[length + 1];
        std::strcpy(seq, fragment.c_str());
        //stopped here.
        translate2numbers((uchar *) seq, (uint32) length, mBwtfmiDB->astruct);
        SI *si = NULL;
        if (num_mm > 0) {
            if (num_mm == mOptions->mDNASearchOptions->misMatches) { //after last mm has been done, we need to have at least reached the min_length
                si = maxMatches_withStart(mBwtfmiDB->fmi, seq, (uint32) length, mOptions->mDNASearchOptions->minFragLength, 1, t->si0, t->si1, t->matchlen);
            } else {
                si = maxMatches_withStart(mBwtfmiDB->fmi, seq, (uint32) length, t->matchlen, 1, t->si0, t->si1, t->matchlen);
            }
        } else {
            si = maxMatches(mBwtfmiDB->fmi, seq, (uint32) length, mOptions->mDNASearchOptions->seedLength, 0); //initial matches
        }

        if (!si) { // no match for this fragment
            if (mOptions->debug) {
                std::stringstream ss;
                ss << "No match for this fragment (DNA).";
                cCout(ss.str(), 'y');
            }
            delete[] seq;
            delete t;
            continue; // continue with the next fragment
        }

        if (mOptions->debug) {
            std::stringstream ss;
            ss << "Longest match has length (DNA) " << (uint) si->ql;
            cCout(ss.str(), 'y');
        }

        if (mOptions->mDNASearchOptions->misMatches > 0 && num_mm < mOptions->mDNASearchOptions->misMatches) {
            SI *si_it = si;
            while (si_it) {
                uint32 match_right_end = si_it->qi + si_it->ql - 1;
                if (num_mm > 0)
                    assert(match_right_end == length - 1); // greedy matches end always at the end
                if (mOptions->debug) {
                    std::stringstream ss;
                    ss << "DNA Match from " << si_it->qi << " to " << match_right_end << ": " << fragment.substr(si_it->qi, match_right_end - si_it->qi + 1) << " (" << si_it->ql << ")";
                    cCout(ss.str(), 'y');
                }
                if (si_it->qi > 0 && match_right_end + 1 >= mOptions->mDNASearchOptions->minFragLength) {
                    //1. match must end before beginning of fragment, i.e. it is extendable
                    //2. remaining fragment, from zero to end of current match, must be longer than minimum length of accepted matches
                    const size_t erase_pos = (match_right_end < length - 1) ? match_right_end + 1 : std::string::npos;
                    addAllMismatchVariantsAtPosSI(t, (uint) (si_it->qi - 1), erase_pos, si_it);
                }
                si_it = si_it->samelen ? si_it->samelen : si_it->next;
            }
        }


        if ((unsigned int) si->ql < mOptions->mDNASearchOptions->minFragLength) { // match was too short
            if (mOptions->debug) {
                std::stringstream ss;
                ss << "Match of length (DNA) " << si->ql << " is too short";
                cCout(ss.str(), 'y');
            }
            delete[] seq;
            delete t;
            recursive_free_SI(si);
            continue; // continue with the next fragment
        }

        eval_match_scores(si, t);

        delete[] seq;
        delete t;
    }

    if (best_matches_SI.empty()) {
        return;
    }

    if (mOptions->mDNASearchOptions->useEvalue) {
        //calc e-value and only return match if > cutoff

        double bitscore = (dLAMBDA * best_match_score - dLN_K) / dLN_2;
        double Evalue = mBwtfmiDB->db_length * query_len * pow(2, -1 * bitscore);
        if (mOptions->debug) {
            std::stringstream ss;
            ss << "E-value = " << Evalue;
            cCout(ss.str(), 'y');
        }

        if (Evalue > mOptions->mDNASearchOptions->minEvalue) {
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
         //recursive_free_SI(itm); //why not this one;
        free(itm);
    }
}

void DNASearcher::eval_match_scores(SI *si, Fragment *frag) {

    if (!si) return;

    // eval the remaining same-length and shorter matches
    if (si->samelen) eval_match_scores(si->samelen, frag);

    if (si->next && si->next->ql >= (int) mOptions->mDNASearchOptions->minFragLength) {
        eval_match_scores(si->next, frag);
    } else if (si->next) {
        recursive_free_SI(si->next);
    }

    unsigned int score = calcScore(frag->seq, si->qi, si->ql, frag->diff);

    if (mOptions->debug) {
        std::stringstream ss;
        ss << "Match " << frag->seq.substr(si->qi, si->ql) << " (length(DNA)=" << (unsigned int) si->ql << " score=" << score << " num_mm=" << frag->num_mm << ")";
        cCout(ss.str(), 'y');
    }
    if (score < mOptions->mDNASearchOptions->minScore) {
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
        if (mOptions->verbose) {
            best_matches.clear();
            best_matches.push_back(frag->seq.substr(si->qi, si->ql));
        }
    } else if (score == best_match_score && best_matches_SI.size() < mOptions->mDNASearchOptions->max_matches_SI) {
        best_matches_SI.push_back(si);
        if (mOptions->verbose)
            best_matches.push_back(frag->seq.substr(si->qi, si->ql));
    } else {
        free(si);
        si = NULL;
    }
}

void DNASearcher::addAllMismatchVariantsAtPosSI(const Fragment *f, unsigned int pos, size_t erase_pos = std::string::npos, SI *si = NULL) {
    assert(mOptions->mDNASearchOptions.mode == dGREEDY);
    assert(pos < erase_pos);
    assert(f->num_mm == 0 || pos < f->pos_lastmm);

    std::string fragment = f->seq; // make a copy to modify the sequence at pos
    assert(fragment.length() >= mOptions->mDNASearchOptions->minFragLength);
    char origchar = fragment[pos];
    assert(blosum_subst.count(origchar) > 0);
    if(blosum_subst.find(origchar) == blosum_subst.end()){
        if(mOptions->verbose){
            std::cerr << "Warning: No substitutions found for nucleotide " << origchar << " in pos: " << pos << " in frag: " << fragment << std::endl;
        }
        return;
    }

    if (erase_pos != std::string::npos && erase_pos < fragment.length()) {
        if (mOptions->debug) {
            std::stringstream ss;
            ss << "Deleting from position (DNA) " << erase_pos;
            cCout(ss.str(), 'y');
        }
        fragment.erase(erase_pos);
    }

    //calc score for whole sequence, so we can substract the diff for each substitution
    unsigned int score = calcScore(fragment, f->diff) - blosum62diag[nuc2int[(uint8_t) origchar]];
    IndexType siarray[2], siarrayupd[2];
    siarray[0] = si->start;
    siarray[1] = si->start + (IndexType) si->len;
    for (auto itv : blosum_subst.at(origchar)) {
        // we know the difference between score of original aa and substitution score, this
        // has to be subtracted when summing over all positions later
        // so we add this difference to the fragment
        int score_after_subst = score + b62[nuc2int[(uint8_t) origchar]][nuc2int[(uint8_t) itv]];
        if (score_after_subst >= (int) best_match_score && score_after_subst >= (int) mOptions->mDNASearchOptions->minScore) {
            if (UpdateSI(mBwtfmiDB->fmi, mBwtfmiDB->astruct->trans[(size_t) itv], siarray, siarrayupd) != 0) {
                fragment[pos] = itv;
                int diff = b62[nuc2int[(uint8_t) origchar]][nuc2int[(uint8_t) itv]] - blosum62diag[nuc2int[(uint8_t) itv]];
                if (mOptions->debug) {
                    std::stringstream ss;
                    ss << "Adding fragment   " << fragment << " with mismatch at pos " << pos << " ,diff " << f->diff + diff << ", max score " << score_after_subst;
                    cCout(ss.str(), 'y');
                }
                fragments.emplace(score_after_subst, new Fragment(fragment, f->num_mm + 1, pos, f->diff + diff, siarrayupd[0], siarrayupd[1], si->ql + 1));
            } else if (mOptions->debug) {
                fragment[pos] = itv;
                std::stringstream ss;
                ss << "Skipping fragment " << fragment << " mismatch at pos " << pos << ", because " << itv << " is not a valid extension";
                cCout(ss.str(), 'y');
            }
        } else {
            if (mOptions->debug) {
                fragment[pos] = itv;
                std::stringstream ss;
                ss << "Skipping fragment " << fragment << " and following fragments, because score is too low: " << score_after_subst << " < " << std::max(best_match_score, mOptions->mDNASearchOptions->minScore);
                cCout(ss.str(), 'y');
            }
            break;
        }
    }
}

Fragment* DNASearcher::getNextFragment(uint min_score) {
    if (fragments.empty()) {
        return NULL;
    }

    auto it = fragments.begin();

    //if (mOptions->debug) std::cerr << "max fragment score/length (DNA) = " << it->first << "\n";
    if (it->first < min_score) { //the highest scoring fragment in the sorted list is below threshold, then search stops
        return NULL;
    }

    Fragment *f = it->second;
    if (mOptions->debug) {
        std::stringstream ss;
        ss << "Fragment (DNA) = " << f->seq << "; max fragment score/length = " << it->first;
        cCout(ss.str(), 'y');
    }
    fragments.erase(it);

    return f;
}

void DNASearcher::ids_from_SI(SI *si) {
    IndexType k, pos;
    int iseq;
    for (k = si->start; k < si->start + si->len; ++k) {
        if (match_ids.size() > mOptions->mDNASearchOptions->max_match_ids) {
            break;
        }
        get_suffix(mBwtfmiDB->fmi, mBwtfmiDB->bwt->s, k, &iseq, &pos);
        match_ids.insert(mBwtfmiDB->bwt->s->ids[iseq]);
    }
}

void DNASearcher::ids_from_SI_recursive(SI *si) {
    SI *si_it = si;
    //if (mOptions->debug) std::cerr << "debug: ids_from_SI_recursive: a" << "\n";
    while (si_it) {
        IndexType k, pos;
        int iseq;
        // if (mOptions->debug) std::cerr << "debug: ids_from_SI_recursive: b" << "\n";
        for (k = si_it->start; k < si_it->start + si_it->len; ++k) {
            // if (mOptions->debug) std::cerr << "debug: ids_from_SI_recursive: c" << "\n";
            //std::cerr << "debug: DNASearcher::ids_from_SI_recursive: c; k: " << k << "; iseq: " << iseq << "; pos: " << pos << "\n";
            if (match_ids.size() > mOptions->mDNASearchOptions->max_match_ids) {
                break;
            }
            //std::cerr << "debug: DNASearcher::ids_from_SI_recursive: d; k: " << k << "; iseq: " << iseq << "; pos: " << pos << "\n";
            get_suffix(mBwtfmiDB->fmi, mBwtfmiDB->bwt->s, k, &iseq, &pos);
            // std::cerr << "debug: DNASearcher::ids_from_SI_recursive: e; k: " << k << "; iseq: " << iseq << "; pos: " << pos << "\n";
            match_ids.insert(mBwtfmiDB->bwt->s->ids[iseq]);
            //std::cerr << "debug: DNASearcher::ids_from_SI_recursive: f; k: " << k << "; iseq: " << iseq << "; pos: " << pos << "\n";
        } // end for
        //std::cerr << "debug: DNASearcher::ids_from_SI_recursive: g" << "\n";
        si_it = si_it->samelen;
        //std::cerr << "debug: DNASearcher::ids_from_SI_recursive: h" << "\n";
    } // end while all SI with same length
    // std::cerr << "debug: DNASearcher::ids_from_SI_recursive: i" << "\n";
}

uint DNASearcher::calcScore(const std::string &s, size_t start, size_t len, int diff) {
    int score = 0;
    for (size_t i = start; i < start + len; ++i) {
        score += blosum62diag[nuc2int[(uint8_t) s[i]]];
    }
    score += diff;
    return score > 0 ? score : 0;
}

uint DNASearcher::calcScore(const std::string &s, int diff) {
    int score = 0;
    for (size_t i = 0; i < s.length(); ++i) {
        score += blosum62diag[nuc2int[(uint8_t) s[i]]];
    }
    score += diff;
    return score > 0 ? score : 0;
}

uint DNASearcher::calcScore(const std::string &s) {
    unsigned int score = 0;
    for (size_t i = 0; i < s.length(); ++i) {
        score += blosum62diag[nuc2int[(uint8_t) s[i]]];
    }
    return score;
}
