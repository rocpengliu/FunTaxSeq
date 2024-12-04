/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */


/* Since the definition of ex2 may differ between implementations,
   this header file MUST be included AFTER defining ex2
*/


#ifndef COMMONFMI_h
#define COMMONFMI_h

// const char ex1=16;   // Exponent for checkpoint index 1 (dist 2^e1 between them)
#define ex1 16   // Exponent for checkpoint index 1 (dist 2^e1 between them)

#define CHUNK_SIZE 1024 * 1024 * 512

/* If ex2=3 the binary numbers are:
   size2 =   ...001000   (2^ex2=8)
   check2 =  ...000111 = (1<<ex2)-1
   round2 = ~...000111 = ...111000 (used to round -> set the lowest ex2 bits to zero)
*/
const IndexType size1 = (IndexType)1 << ex1;
//const IndexType check1 = size1-1;
const IndexType check1 = ((IndexType)1 << ex1)-1;
const IndexType size2 = (IndexType)1 << ex2;
// const IndexType check2 = size2-1;
const IndexType check2 = ((IndexType)1 << ex2)-1;
// const IndexType round2 = ~(size2-1);          /* 1 at every bit above bit ex2-1 */
const IndexType round2 = ~(((IndexType)1 << ex2)-1);          /* 1 at every bit above bit ex2-1 */



/*
  Find length to next index2 or end of BWT from position k.
*/
static inline int fmi_end_length(IndexType k, IndexType bwtlen) {
  IndexType k2 = k&round2;
  if ( k2+size2>bwtlen ) return bwtlen-k;
  else return k2+size2-k;
}



/* Where is the nearest chkpt2? Forward (+1) or backward (-1)
   If the bit ex2-1 is set, the number is half way or more from the previous
   checkpoint. So with mask = 1 << (ex2-1) the check is k&mask
*/
static inline int fmi_direction(IndexType k) {
  const IndexType mask = (IndexType)1 << (ex2-1);
  if ( k&mask ) return 1;
  else return -1;
}



/* Get the checkpointed FMI value for k and direction
   If direction=-1, the checkpoints are k>>ex2 and k>>ex1
   if          = 1, the checkpoints are chpt2=k>>ex2+1 and chpt1 = chpt2>>(ex1-ex2)
*/
static inline IndexType fmi_chpt_value_with_dir(const FMI *f, const IndexType k, const uchar c, const int direction) {
  IndexType chpt1, chpt2;

  chpt2 = k>>ex2;
  if (direction>0) chpt2 += 1;

  /* Note that if direction is +1 and we're just below a checkpoint1, then
     the checkpoint 1 to use is (k>>ex1)+1 rather than k>>ex1
     We can always use chpt1=chpt2>>(ex1-ex2)
  */
  chpt1 = chpt2>>(ex1-ex2);

  return f->index1[chpt1][c] + f->index2[chpt2][c];
}





/* Methods differ in allocation for index2
   + some additionals
*/
static FMI *alloc_FMI_common(uchar *bwt, IndexType bwtlen, int alen, size_t index2_size) {
  int i;
  FMI *f = (FMI*)malloc(sizeof(FMI));
  f->alen = alen;
  f->bwt = bwt;
  f->bwtlen = bwtlen;
  f->N1 = ((bwtlen-1)>>ex1)+2;
  if (f->N1<<ex1 == bwtlen) f->N1 -= 1;
  f->N2 = ((bwtlen-1)>>ex2)+2;
  if (f->N1<<ex2 == bwtlen) f->N2 -= 1;
  f->index1 = (IndexType **)malloc(f->N1*sizeof(IndexType *));
  for (i=0;i<f->N1;++i) f->index1[i]=(IndexType *)malloc(alen*sizeof(IndexType));
  f->index2 = (ushort**)malloc(f->N2*sizeof(ushort*));
  for (i=0;i<f->N2;++i) f->index2[i]=(ushort *)malloc(index2_size*alen);
  return f;
}



/* This function sets the values at index1 and index2.
   Each FMI method may have to additionally recode the BWT
*/
static FMI *makeIndex_common(uchar *bwt, long bwtlen, int alen) {
  IndexType i, ii, R1, R2, *total;
  int a, *current;
  // uchar *sbwt;
  FMI *fmi = alloc_FMI(bwt,bwtlen,alen);

  current = (int *)calloc(alen,sizeof(int));
  total = (IndexType *)calloc(alen,sizeof(IndexType));
  for (a=0;a<alen;++a) fmi->index2[0][a]=0;

  /*
  for (a=0;a<alen;++a) fmi->index1[0][a]=0;
  for (a=0;a<alen;++a) current[a]=total[a]=0;
  */

  // Note that current char is NOT counted
  i=0;
  R1=R2=0;
  // sbwt=bwt;
  IndexType ten_percent = bwtlen/10;
  for (ii=0; ii<bwtlen; ++ii) {
    if ( (ii+1)%ten_percent==0 ) fprintf(stderr,"%d%% ... ",10*(int)((ii+1)/ten_percent));
    /* Check if we are at a checkpoint 1 */
    if ( !(ii&check1) ) {
      R1 = ii>>ex1;
      for (a=0;a<alen;++a) fmi->index1[R1][a]=total[a];
    }
    /* Check if we are at a checkpoint 2 */
    if ( ii>0 && !(ii&check2) ) {
      R2 = ii>>ex2;
      /* Checkpoint values */
      for (a=0; a<alen; ++a) fmi->index2[R2][a]=(ushort)(total[a]-fmi->index1[R1][a]);
      /* Reset counter */
      for (a=0;a<alen;++a) current[a]=0;
      i=0;
      // sbwt+=size2;
    }
    // a = sbwt[i];
    a = bwt[ii];
    if (a<0 || a>=alen) {
      fprintf(stderr,"makeIndex_common: letter %d not in range at %ld, alen=%d\n",a,ii,alen);
      exit(199);
    }
    total[a] += 1;
    current[a] += 1;
    // ++i;
  }

  // The last entry of index2 (doesn't matter if last was already a power of ex2)
  // R2 = 1+((ii-1)>>ex2);
  // R2 = 1+(ii>>ex2);
  R2 = fmi->N2-1;
  /* Insert values in checkpoint 2 */
  for (a=0;a<alen;++a) fmi->index2[R2][a]=(ushort)(total[a]-fmi->index1[R1][a]);

  fprintf(stderr,"index2 done ... ");

  // Save the letter starts in the last index1
  // Add to all of index1
  fmi->index1[fmi->N1-1][0]=0;
  for (a=1;a<alen;++a) fmi->index1[fmi->N1-1][a]=fmi->index1[fmi->N1-1][a-1]+total[a-1];
  for (R1=0; R1 < fmi->N1-1; ++R1) for (a=1;a<alen;++a) fmi->index1[R1][a] += fmi->index1[fmi->N1-1][a];

  free(current);
  free(total);

  return fmi;
}


/* Write the FMI in file (binary) */
static void write_fmi_common(const FMI *f, int index2_size, FILE *fp) {
  int i;
  fwrite(&(f->alen),sizeof(int),1,fp);
  fwrite(&(f->bwtlen),sizeof(IndexType),1,fp);
  fwrite(&(f->N1),sizeof(int),1,fp);
  fwrite(&(f->N2),sizeof(int),1,fp);
  fwrite(f->bwt,sizeof(uchar),f->bwtlen,fp);
  for (i=0; i<f->N1; ++i) fwrite(f->index1[i],sizeof(IndexType),f->alen,fp);
  for (i=0; i<f->N2; ++i) fwrite(f->index2[i],1,f->alen*index2_size,fp);
}



/* Read the FMI in file (binary)
*/
static FMI *read_fmi_common(int index2_size, FILE *fp) {
  printf("starting to Reading FMI common\n");
  int i;
  FMI *f = (FMI *)malloc(sizeof(FMI));
  printf("allocating memory for reading FMI complete\n");
  f->bwt=NULL;

  printf("reading alen, bwtlen, N1, N2 complete\n");
  fread(&(f->alen),sizeof(int),1,fp);
  fread(&(f->bwtlen),sizeof(IndexType),1,fp);
  fread(&(f->N1),sizeof(int),1,fp);
  fread(&(f->N2),sizeof(int),1,fp);
  printf("bwt and index1 reading complete\n");

  printf("allocating memory for bwt and index1\n");
  f->bwt=(uchar *)malloc(f->bwtlen*sizeof(uchar));
  printf("allocating memory for bwt and index1 completed\n");
  printf("reading bwt len %ld\n", f->bwtlen);

  if(f->bwt==NULL){
    fprintf(stderr,"Error: malloc failed for f->bwt\n");
    exit(1);
  }

  //######################
  printf("starting to read bwt by chuncks\n");
  size_t bytes_read = 0;
  while(bytes_read < f->bwtlen){
    size_t chunk_size = ((f->bwtlen - bytes_read) > CHUNK_SIZE) ? CHUNK_SIZE : (f->bwtlen - bytes_read);
    printf("Reading chunk of size %zu, remaining bytes: %zu\n", chunk_size, f->bwtlen - bytes_read);
    size_t n = fread(f->bwt + bytes_read, sizeof(uchar), chunk_size, fp);
    if (n < chunk_size) {
        if (n == 0) {
            fprintf(stderr, "Error: fread failed for bwt (unexpected EOF or read error)\n");
        } else {
            fprintf(stderr, "Error: fread read %zu bytes, expected %zu\n", n, chunk_size);
        }
        exit(1);
    }
    bytes_read += n;
    printf("Bytes read so far: %zu\n", bytes_read);
  }
//######################
  //fread(f->bwt, sizeof(uchar), f->bwtlen, fp); // why can not read and killed. this is the oringal
  printf("f->bwt reading completed with reading length %zu\n", bytes_read);

  f->index1 = (IndexType **)malloc(f->N1*sizeof(IndexType *));
  printf( "f->index1 allocation complete\n");
  for (i=0;i<f->N1;++i) {
    f->index1[i]=(IndexType *)malloc(f->alen*sizeof(IndexType));
    fread(f->index1[i],sizeof(IndexType),f->alen,fp);
  }
  printf("index1 reading complete\n");




// Allocate memory for all rows in a single block
  size_t total_rows_size = (size_t)f->N2 * (size_t)f->alen * (size_t)index2_size;
  f->index2 = (ushort **)malloc(f->N2 * sizeof(ushort *));
  ushort *block = (ushort *)malloc(total_rows_size);

  if (f->index2 == NULL || block == NULL) {
      fprintf(stderr, "Memory allocation failed for index2\n");
      exit(1);
  }

  // Assign pointers to each row
  for (i = 0; i < f->N2; ++i) {
      f->index2[i] = block + (i * f->alen * index2_size / sizeof(ushort));
  }

  // Read the entire block
  size_t bytes_read = fread(block, 1, total_rows_size, fp);
  if (bytes_read < total_rows_size) {
      fprintf(stderr, "Error: fread failed for index2 block (read %zu bytes instead of %zu)\n",
              bytes_read, total_rows_size);
      free(block);
      free(f->index2);
      exit(1);
  }



//   //****************************************** */
//   printf("allocating memory for N2 with N2 is %d\n", f->N2);
//   f->index2 = (ushort **)malloc(f->N2 * sizeof(uchar *));//orignal code
//   if (f->index2 == NULL) {
//     fprintf(stderr, "Memory allocation for index2 failed\n");
//     exit(1);
//   }
//   printf("allocating memory for N2 completed\n");

//   size_t total_memory = (size_t)f->N2 * (size_t)f->alen * (size_t)index2_size * sizeof(ushort);
//   printf("Total memory for index2: %zu bytes\n", total_memory);

//   for (i=0;i<f->N2;++i) {
//     f->index2[i] = (ushort *)malloc(f->alen * index2_size);
//      if (f->index2[i] == NULL) {
//         fprintf(stderr, "Memory allocation for index2[%d] failed\n", i);
//         exit(1);
//     }
//     //fread(f->index2[i],1,f->alen*index2_size,fp);//original

// //######################
//     size_t bytes_to_read = f->alen * index2_size;
//     size_t bytes_read = fread(f->index2[i], 1, bytes_to_read, fp);
//     if (bytes_read < bytes_to_read) {
//         fprintf(stderr, "Error: fread failed for f->index2[%d] (read %zu bytes instead of %zu)\n",
//                     i, bytes_read, bytes_to_read);
//         exit(1);
//     }

// //######################

//     if (i % 10000000 == 0) {
//         printf("Reading index2 row %d\n", i);
//     }
//   }


  //******************************************* */
  printf("index2 reading complete\n");
  return f;
}

#endif
