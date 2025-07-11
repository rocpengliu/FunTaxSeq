/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "sequence.h"


SEQstruct *alloc_SEQstruct() {
  SEQstruct *ss=(SEQstruct *)malloc(sizeof(SEQstruct));
  ss->len = 0;
  ss->rc = 0;
  ss->pos = 0;
  ss->id = NULL;
  ss->descr = NULL;
  ss->start = NULL;
  ss->id_filepos = 0;
  ss->seq_filepos = 0;
  ss->sort_order = 0;
  ss->next=NULL;
  return ss;
}

void free_SEQstruct(SEQstruct *ss) {
  if (ss) {
    if (ss->id) {
      free(ss->id);
      ss->id = NULL;
    }
    if(ss->descr) {
      free(ss->descr);
      ss->descr = NULL; 
    }            // added;
    if(ss->start) {
      free(ss->start);//added;
      ss->start = NULL;
    }
    free(ss);
    ss = NULL;
  }
}

/* Assumes that base->start points to the whole sequence */
void recursive_free_SEQstruct(SEQstruct *base) {
  SEQstruct *ss, *next;
  if (base->start){//redundant, because free_SEQstruct will free start
      free(base->start);
      base->start = NULL;
  }

  ss=base;
  next=ss->next;
  while (ss) {
    free_SEQstruct(ss);
    ss=next;
    if (next) next=ss->next;
  }
}

/* Makes a translation table from an alphabet to a translation, so
        table[alphabet[i]] = translation[i]

   If translation==NULL, the letter alphabet[i] is translated to i.

   Letters not in alphabet are translated to dummy.
   If dummy==0, dummy is set to the translation of the last char in alphabet
   (assumed to be a "wildcard" character)

   Translation for non-characters is -1.

   casesens !=0, means case sensitive, otherwise case insentitive.

   Returns an array (char *table) of length 128

   The length of the translation table (which may contain zeros)
   HAS to be as long (or longer) than alphabet.

   If 0 is not in the alphabet, it is translated to 0

   peng added
     this is to convert the letter alphabet to number, it liletrly convert the ASCII table 
  to the corresponding number from 0:len(alphabet), say 0:3 and 0:19 for DNA and aa, the others will be -1
  for no alphabet and 4 or 21 for alphabet which are not in the DNA or aa alphabet.
  eg
  out[i0] = 0
  out[i1] = -1
  out[i41] = -1
  out[i42] = 4
  out[i43] = -1
  out[i64] = -1
  out[i65] = 0 A
  out[i66] = 4
  out[i67] = 1 C
  out[i68] = 4
  out[i71] = 2 G
  out[i72] = 4
  out[i84] = 3 T
  out[i85] = 4
  out[i91] = -1
  out[i97] = 0
  out[i98] = 4
  out[i99] = 1
  out[i100] = 4
  out[i103] = 2
  out[i104] = 4
  out[i116] = 3
  out[i117] = 4
  out[i123] = -1
*/
static char *translation_table(char *alphabet, char *translation, char dummy, int casesens) {
  int i, l, freetrans=0;
  char *table = (char*)malloc(128*sizeof(char));//128

  l = strlen(alphabet);
  //printf("translation_table alphabet alen: %d\n", l);
  if (translation==0) {
    translation = (char *)malloc(l*sizeof(char));
    for (i=0; i<l; ++i) translation[i]=i;
    freetrans=1;
  }
  if (dummy==0) dummy = translation[l-1];

  table[0]=0;
  for (i=1; i<128; ++i) {
    table[i]=-1;
    if (isalpha((char)i)) table[i]=dummy;
  }

  if (!casesens) {
    for (i=0; i<l; ++i) {
      table[toupper(alphabet[i])]=translation[i];
      table[tolower(alphabet[i])]=translation[i];
    }
  } else {
    for (i=0; i<l; ++i) table[alphabet[i]] = translation[i];
  }
  if (freetrans) free(translation);
  return table;
}

static char *dnaComplement(char *alphabet) {
  int l=strlen(alphabet);
  char *comp=(char *)malloc((l+1)*sizeof(char));
  int i;
  for (i=0;i<l;++i) {
    switch (alphabet[i]) {
        case 'a': comp[i]='t'; break;
        case 'A': comp[i]='T'; break;
        case 't': comp[i]='a'; break;
        case 'T': comp[i]='A'; break;
        case 'c': comp[i]='g'; break;
        case 'C': comp[i]='G'; break;
        case 'g': comp[i]='c'; break;
        case 'G': comp[i]='C'; break;
        default: comp[i]=alphabet[i];
    }
  }
  comp[i]=0;
  return comp;
}

AlphabetStruct *alloc_AlphabetStruct(char *a, int caseSens, int revcomp) {
  AlphabetStruct *astruct = (AlphabetStruct *)malloc(sizeof(AlphabetStruct));
  astruct->a = strdup(a);
  astruct->len = strlen(a);
  astruct->caseSens = caseSens;
  astruct->trans = translation_table(a, NULL, astruct->len-1, astruct->caseSens);
  astruct->comp = (revcomp ? dnaComplement(a) : NULL);
  return astruct;
}

void free_AlphabetStruct(AlphabetStruct *astruct) {
  if (astruct) {
    if (astruct->a){
      free(astruct->a); astruct->a = NULL;
    }
    if (astruct->trans){
      free(astruct->trans); astruct->trans = NULL;
    }
    if (astruct->comp) {
      free(astruct->comp);astruct->comp = NULL;
    }
    free(astruct);
    astruct = NULL;
  }
}

/*
  translate a sequence (s) to numbers
 */
void translate2numbers(uchar *s, const IndexType slen, AlphabetStruct *astruct) {
  IndexType k;
  for (k=0;k<slen;++k){
      s[k]=astruct->trans[s[k]];
  }
}