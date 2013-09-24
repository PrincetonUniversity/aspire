// Original readdata routines are in orig_readdata.c

#include "myinclude.h"
#ifdef __MEX
#include "mex.h"
#endif

// #define MYCALLOC(VAR,TYPE,SIZE) VAR = (TYPE*)calloc(SIZE, sizeof(TYPE))
// #define MYFREE(VAR) free(VAR)

// #define DATABLOCKIND(DATA,BLOCK,NUMBLOCK) ((DATA+1)-1)*NUMBLOCK + BLOCK - 1

size_t quicksort5(size_t*, size_t*, size_t*, size_t*, double*, size_t, size_t);
size_t partition5(size_t*, size_t*, size_t*, size_t*, double*, size_t, size_t);
void skip_to_end_of_line(FILE*);
size_t get_line(FILE*, char*, size_t);
size_t max_line_length(FILE*);

extern void dsyev_();

size_t readdata_sdpa(char* datafilename, size_t* passed_m, size_t* passed_numblk, ptrdiff_t** passed_blksz,
                  char** passed_blktype, double** passed_b, double** passed_CAent,
                  size_t** passed_CArow, size_t** passed_CAcol, size_t** passed_CAinfo_entptr,
                  size_t** passed_CAinfo_rowcolptr, char** passed_CAinfo_type,
                  char** passed_CAinfo_storage)
{
  // declare variables 
  size_t i, ct, buflen, ret, *nnz, *temp_num, *temp_blk, specindex;
  int num, blk, pos1, pos2;
  char *buf, c, *ptr1, *ptr2;
  double entry;
  FILE *datafile;

  size_t sam_nnz, *sam_num, *sam_blk, *sam_row, *sam_col;
  size_t tmp_nnz, *tmp_num, *tmp_blk, *tmp_row, *tmp_col;
  double *sam_ent;
  double *tmp_ent;
  size_t tmp_ct;
  int flag;

  int m, numblk;
  ptrdiff_t *blksz;
  char   *blktype;
  double *b, *CAent;
  size_t    *CArow, *CAcol, *CAinfo_entptr, *CAinfo_rowcolptr;
  char   *CAinfo_type, *CAinfo_storage;


  ////////////////////////////////
  // Step 0
  ////////////////////////////////

  // Open file, and determine the length of the longest line. 
  // Note: this can fail if last line doesn't have a new line character; how to fix? 

  datafile = fopen(datafilename, "r");
  if(datafile == NULL) {
    printf("Can't get file %s\n", datafilename);
    exit(0);
  }
  buflen = max_line_length(datafile) + 10;
  fclose(datafile);

  // Allocate buffer space, which will hold each line that is read in

  MYCALLOC(buf, char, buflen);


  ////////////////////////////////
  // Step 1
  ////////////////////////////////


  // Get m, numblk, *blksz, *blktype, *b;
  // size of partitions of (CAent, CArow, CAcol) into (data matrix)-(block) pairs

  // Open file

  datafile = fopen(datafilename, "r");

  // Read through the comment lines. 
 
  c = getc(datafile);
  while(c == '"' || c == '*') {
    skip_to_end_of_line(datafile);
    c = getc(datafile);
  }
  ungetc(c, datafile);

  // Get m

  ret = get_line(datafile, buf, buflen);
  sscanf(buf, "%d", &m);

  // Get numblk

  ret = get_line(datafile, buf, buflen);
  sscanf(buf, "%d", &numblk);

  // Prepare to get blksz and blktype
 
  MYCALLOC(blksz, ptrdiff_t, numblk);
  MYCALLOC(blktype, char, numblk);

  // Get blksz

  ret = get_line(datafile, buf, buflen);
  ptr1 = buf;
	for(i = 0; i < numblk; i++) {
		blksz[i] = strtol(ptr1, &ptr2, 10);
	  ptr1 = ptr2;
  }

  // Get blktype (easy from blksz)
  // Note: In SDPA format, X nonneg variables constitutes
  // one (diagonal) block of size X, which is indicated by
  // a negative block size.

  for(i = 0; i < numblk; i++) {
    if(blksz[i] < 0) {
      blktype[i] = 'd';
      blksz[i] = -blksz[i];
    }
    else blktype[i] = 's';
  }

  // Prepare to get b

  MYCALLOC(b, double, m);

  // Get b 

  ret = get_line(datafile, buf, buflen);
  ptr1 = buf;
  for(i = 0; i < m; i++) {
	  b[i] = strtod(ptr1, &ptr2);
	  ptr1 = ptr2;
	}

  // Count number of remaining lines

  sam_nnz = 0;
  ret = fscanf(datafile, "%d %d %d %d %le", &num, &blk, &pos1, &pos2, &entry);
  do {
    sam_nnz++;
    ret = fscanf(datafile, "%d %d %d %d %le", &num, &blk, &pos1, &pos2, &entry);
  } while(ret == 5);
  fclose(datafile);

  // Store remaining lines (re-reading file is similar to above)
  
  MYCALLOC(sam_num, size_t, sam_nnz+1); // Note: length+1 vectors
  MYCALLOC(sam_blk, size_t, sam_nnz+1);
  MYCALLOC(sam_row, size_t, sam_nnz+1);
  MYCALLOC(sam_col, size_t, sam_nnz+1);
  MYCALLOC(sam_ent, double, sam_nnz+1);

  datafile = fopen(datafilename, "r");
  c = getc(datafile);
  while(c == '"' || c == '*') {
    skip_to_end_of_line(datafile);
    c = getc(datafile);
  }
  ungetc(c, datafile);
  ret = get_line(datafile, buf, buflen);
  ret = get_line(datafile, buf, buflen);
  ret = get_line(datafile, buf, buflen);
  ret = get_line(datafile, buf, buflen);

  flag = 0;
  for(ct = 1; ct <= sam_nnz; ct++) {
    ret = fscanf(datafile, "%d %d %d %d %le", &num, &blk, &pos1, &pos2, &entry);
    if(pos1 == 0 || pos2 == 0) {
      printf("Error (readdata_sdpa): Encountered '0' row or column index.\n");
      exit(0);
    }
    sam_num[ct] = num;
    sam_blk[ct] = blk;
    if(pos1 <= pos2) {
      sam_row[ct] = pos1;
      sam_col[ct] = pos2;
    }
    else { 
      flag = 1;
      sam_row[ct] = pos2;
      sam_col[ct] = pos1;
    }
    sam_ent[ct] = entry;
  }
  fclose(datafile);
  if(flag) printf("Warning: Switching lower-triangular entries to upper-triangular.\n");

  quicksort5(sam_num, sam_blk, sam_row, sam_col, sam_ent, 0, sam_nnz);

  // Loop through entries and accumulate repeated entries.
  // Some entries become zero.
  // Is this correct? Think so.

  flag = 0;
  do {
    tmp_ct = 0;
    for(ct = 2; ct <= sam_nnz; ct++)
      if(sam_num[ct] == sam_num[ct-1] &&
         sam_blk[ct] == sam_blk[ct-1] &&
         sam_row[ct] == sam_row[ct-1] &&
         sam_col[ct] == sam_col[ct-1] &&
         sam_ent[ct] != 0.0) {
        sam_ent[ct-1] += sam_ent[ct];
        sam_ent[ct] = 0.0;
        flag = 1;
        tmp_ct++;
      }
  } while(tmp_ct > 0);
  if(flag) printf("Warning: Accumulating value of repeated upper-triangular entry.\n");

  // Delete zero entries, some of which may have been created above.
  // Requires temp space.

  tmp_nnz = 0;
  for(ct = 1; ct <= sam_nnz; ct++)
    if(fabs(sam_ent[ct]) > DBL_EPSILON)
      tmp_nnz++;

  MYCALLOC(tmp_num, size_t, tmp_nnz+1);
  MYCALLOC(tmp_blk, size_t, tmp_nnz+1);
  MYCALLOC(tmp_row, size_t, tmp_nnz+1);
  MYCALLOC(tmp_col, size_t, tmp_nnz+1);
  MYCALLOC(tmp_ent, double, tmp_nnz+1);

  tmp_ct = 1;
  for(ct = 1; ct <= sam_nnz; ct++) {
    if(fabs(sam_ent[ct]) > DBL_EPSILON) {
      tmp_num[tmp_ct] = sam_num[ct];
      tmp_blk[tmp_ct] = sam_blk[ct];
      tmp_row[tmp_ct] = sam_row[ct];
      tmp_col[tmp_ct] = sam_col[ct];
      tmp_ent[tmp_ct] = sam_ent[ct];
      tmp_ct++;
    }
  }

  MYFREE(sam_num);
  MYFREE(sam_blk);
  MYFREE(sam_row);
  MYFREE(sam_col);
  MYFREE(sam_ent);

  sam_num = tmp_num;
  sam_blk = tmp_blk;
  sam_row = tmp_row;
  sam_col = tmp_col;
  sam_ent = tmp_ent;
  sam_nnz = tmp_nnz;

  quicksort5(sam_num, sam_blk, sam_row, sam_col, sam_ent, 0, sam_nnz);

  // Prepare for next step
  
  MYCALLOC(nnz, size_t, (m+1)*numblk);

  // Determine how many nnz entries are in each (data matrix)-(block) pair

  for(ct = 1; ct <= sam_nnz; ct++) {
    if(fabs(sam_ent[ct]) > DBL_EPSILON) nnz[ DATABLOCKIND(sam_num[ct],sam_blk[ct],numblk) ]++;
  }

  ////////////////////////////////
  // Step 2
  ////////////////////////////////

  // Use information gathered from Step 1 to allocate space
  // for remaining structures and to set CAinfo...

  // Allocate all CAinfo...

  MYCALLOC(CAinfo_entptr,    size_t,  (m+1)*numblk + 1);
  MYCALLOC(CAinfo_rowcolptr, size_t,  (m+1)*numblk + 1);
  MYCALLOC(CAinfo_type,      char, (m+1)*numblk);
  MYCALLOC(CAinfo_storage,   char, (m+1)*numblk);

  // Setup all CAinfo...

  ct = 0;

  for(i = 0; i <= m; i++)
    for(blk = 1; blk <= numblk; blk++) {
      
      specindex = DATABLOCKIND(i,blk,numblk);
      CAinfo_entptr[specindex] = ct;
      CAinfo_rowcolptr[specindex] = ct;
      ct += nnz[specindex];

      CAinfo_storage[specindex] = 's'; // SDPA format implies same for all
      if(blktype[blk-1] == 'd') CAinfo_type[specindex] = 'd';
      else CAinfo_type[specindex] = 's';

    }

  CAinfo_entptr[(m+1)*numblk] = ct;
  CAinfo_rowcolptr[(m+1)*numblk] = ct;

  // Allocate remaining space plus some temp space for use in Step 3

  MYCALLOC(CAent, double, ct);
  MYCALLOC(CArow, size_t, ct);
  MYCALLOC(CAcol, size_t, ct);
  MYCALLOC(temp_num, size_t, ct);
  MYCALLOC(temp_blk, size_t, ct);

  ////////////////////////////////
  // Step 3
  ////////////////////////////////

  // Alter for objective
  
  for(ct = 1; ct <= sam_nnz; ct++) {
    temp_num[ct-1] = sam_num[ct];
    temp_blk[ct-1] = sam_blk[ct];
    CArow[ct-1] = sam_row[ct];
    CAcol[ct-1] = sam_col[ct];
    if(sam_num[ct] == 0) CAent[ct-1] = -sam_ent[ct];
    else                 CAent[ct-1] =  sam_ent[ct];
  }

  MYFREE(sam_num);
  MYFREE(sam_blk);
  MYFREE(sam_row);
  MYFREE(sam_col);
  MYFREE(sam_ent);

  ////////////////////////////////
  // Step 4
  ////////////////////////////////

  // Sort CAent, CArow, CAcol according to (data matrix) -> (block) order

  // Not totally sure this is correct but will go with it for now
  quicksort5(temp_num, temp_blk, CArow, CAcol, CAent, 0, CAinfo_entptr[(m+1)*numblk]-1);

  /////////////////////////////////
  // End Step 3
  /////////////////////////////////

  /////////////////////////////////
  // End Step 2
  /////////////////////////////////

  MYFREE(temp_num);
  MYFREE(temp_blk);

  /////////////////////////////////
  // End Step 1
  /////////////////////////////////

  MYFREE(nnz);

  /////////////////////////////////
  // End Step 0
  /////////////////////////////////

  MYFREE(buf);

  // Return back to calling function

  *passed_m                = m;
  *passed_numblk           = numblk;
  *passed_blksz            = blksz;
  *passed_blktype          = blktype;
  *passed_b                = b;
  *passed_CAent            = CAent;
  *passed_CArow            = CArow;
  *passed_CAcol            = CAcol;
  *passed_CAinfo_entptr    = CAinfo_entptr;
  *passed_CAinfo_rowcolptr = CAinfo_rowcolptr;
  *passed_CAinfo_type      = CAinfo_type;
  *passed_CAinfo_storage   = CAinfo_storage;

  return 0;

}

size_t quicksort5(size_t* A1, size_t* A2, size_t* A3, size_t* A4, double* A5, size_t p, size_t r)
{
   size_t q;

   if(p < r) 
   {
      q = partition5(A1, A2, A3, A4, A5, p, r);
      quicksort5(A1, A2, A3, A4, A5, p, q);
      quicksort5(A1, A2, A3, A4, A5, q+1, r);
   }

   return 1;
}

size_t partition5(size_t* A1, size_t* A2, size_t* A3, size_t* A4, double* A5, size_t p, size_t r)
{
   size_t i, j;
   size_t sv1, sv2;
   size_t t1, t2, t3, t4;
   double t5;

   sv1 = A1[p]; sv2 = A2[p];
   i = p-1;
   j = r+1;
   while(i < j) {
      do j--;
      while(A1[j] > sv1 || (A1[j] == sv1 && A2[j] > sv2) );
      do i++;
      while(A1[i] < sv1 || (A1[i] == sv1 && A2[i] < sv2) );
      if(i < j) {
         t1 = A1[j]; t2 = A2[j]; t3 = A3[j]; t4 = A4[j]; t5 = A5[j];
         A1[j] = A1[i]; A2[j] = A2[i]; A3[j] = A3[i]; A4[j] = A4[i]; A5[j] = A5[i];
         A1[i] = t1; A2[i] = t2; A3[i] = t3; A4[i] = t4; A5[i] = t5;
      }
      else return j;
   }

   return 0;
}


// This routine skips to the end of the current line of input from the file datafile. 

void skip_to_end_of_line(FILE *datafile)
{
  char c;
  
  c = getc(datafile);
  while (c != '\n') c = getc(datafile);
}

// This routine reads a line of input size_to a buffer, and translates all
//   occurences of "," "{" "}" "(" ")" to spaces. 

size_t get_line(FILE *datafile, char *buffer, size_t bufsiz)
{
  size_t i, k=0;
  char c;
  
  c = getc(datafile);
  while (c != '\n') {
    buffer[k] = c;
    k++;
    c = getc(datafile);
    if(c == EOF) return(2);
    if(k >= bufsiz) {
      printf("Line too long in input file!  Adjust BUFFERSIZ in readprob.c\n");
      return(1);
    }
  }
  buffer[k] = '\n';
  buffer[k+1] = 0;

  for(i = 0; i <= k; i++) {
    if(buffer[i] == ',' || buffer[i] == '{' ||
       buffer[i] == '}' || buffer[i] == '(' ||
       buffer[i]==')') buffer[i]=' ';
  }

  return(0);
}

size_t max_line_length(FILE *datafile)
{
  size_t maxlen=0, k=0, c;
  
  c = getc(datafile);
  while(c != EOF) {
    k = 0;
    while(c != '\n') {
      c = getc(datafile);
      k++;
    }
    if(k > maxlen) maxlen=k;
    c = getc(datafile);
  }
  
  return(maxlen);
}



size_t readdata_sdplr(char* datafilename, size_t* passed_m, size_t* passed_numblk, ptrdiff_t** passed_blksz,
                   char** passed_blktype, double** passed_b, double** passed_CAent,
                   size_t** passed_CArow, size_t** passed_CAcol, size_t** passed_CAinfo_entptr,
                   size_t** passed_CAinfo_rowcolptr, char** passed_CAinfo_type,
                   char** passed_CAinfo_storage)
{
  // This code assumes a very specific structure of input file,
  // and no error checking is done at all.

  size_t i, j, k, pos1, pos2, num, blk, sz;
  ptrdiff_t *nnz;
  size_t ct, specindex;
  double entry;
  char type;
  FILE *datafile;

  int    m, numblk;
  ptrdiff_t *blksz;
  size_t    ret;
  char   *blktype;
  double *b, *CAent;
  size_t    *CArow, *CAcol, *CAinfo_entptr, *CAinfo_rowcolptr;
  char   *CAinfo_type, *CAinfo_storage;

  ////////////////////////////////
  // Step 1
  ////////////////////////////////

  // Open file
  datafile = fopen(datafilename, "r");
  if(datafile == NULL) {
    printf("Error (readdata_sdplr): Can't get file %s.\n", datafilename);
    exit(0);
  }

  // Read m
  ret = fscanf(datafile, "%d", &m);

  // Read numblk
  ret = fscanf(datafile, "%d", &numblk);

  //if(numblk != 1) { printf("SDPLR data file has more than one block!\n"); exit(0); }

  // Setup blksz and blktype
  MYCALLOC(blksz, ptrdiff_t, numblk);
  MYCALLOC(blktype, char, numblk);

  // Read block size and set blktype
  for(blk = 0; blk < numblk; blk++) {
    ret = fscanf(datafile, "%d", &(blksz[blk]));
    if(blksz[blk] > 0) blktype[blk] = 's';
    else if(blksz[blk] < 0) {
      blksz[blk] *= -1;
      blktype[blk] = 'd';
    }
    else { printf("Problem reading data. Block size 0!\n"); exit(0); }
  }    
  

  // Allocate space for b
  MYCALLOC(b, double, m);

  // Read b
  for(i = 0; i < m; i++)
    ret = fscanf(datafile, "%lf", &(b[i]));

  // Read eval-adjust-val
  // Right now, we will just discard this
  ret = fscanf(datafile, "%lf", &entry);

  // Prepare for next step  
  MYCALLOC(nnz, ptrdiff_t, (m+1)*numblk);

  // Determine how many nnz entries are in each (data matrix)-(block) pair
  // This needs work!
  for(i = 0; i <= m; i++) { 
    for(k = 0; k < numblk; k++) {
      ret = fscanf(datafile, "%d %d %c %d\n", &num, &blk, &type, &sz);
      if(ret < 4) {
        printf("error with fscanf!\n");
        exit(0);
      }
      if(type == 's') {
        nnz[ DATABLOCKIND(num,blk,numblk) ] = 0;
        for(j = 1; j <= sz; j++) {
          ret = fscanf(datafile, "%d %d %lf", &pos1, &pos2, &entry);
          if(fabs(entry) > DBL_EPSILON) nnz[ DATABLOCKIND(num,blk,numblk) ]++;
        }
      }
      else if(type == 'l') {
        nnz[ DATABLOCKIND(num,blk,numblk) ] = -sz*(blksz[blk-1] + 1); // neg num will allow to identify low-rank data matrices later
        for(j = 1; j <= sz*(blksz[blk-1] + 1); j++)
          ret = fscanf(datafile, "%lf", &entry);
          if(ret < 1) {
            printf("error with fscanf!\n");
            exit(0);
          }
      }
    }
  }

  // Finished with this pass through file
  fclose(datafile);

  ////////////////////////////////
  // Step 2
  ////////////////////////////////

  // Allocate all CAinfo...
  MYCALLOC(CAinfo_entptr,    size_t,  (m+1)*numblk + 1);
  MYCALLOC(CAinfo_rowcolptr, size_t,  (m+1)*numblk + 1);
  MYCALLOC(CAinfo_type,      char, (m+1)*numblk);
  MYCALLOC(CAinfo_storage,   char, (m+1)*numblk);

  // Setup all CAinfo...
  ct = 0;

  for(i = 0; i <= m; i++)
    for(blk = 1; blk <= numblk; blk++) {

      specindex = DATABLOCKIND(i,blk,numblk);
      CAinfo_entptr[specindex] = ct;
      CAinfo_rowcolptr[specindex] = ct;
      if(nnz[specindex] > 0) ct += nnz[specindex]; // sparse matrix
      else ct += -nnz[specindex]; // low-rank matrix

      if(nnz[specindex] > 0) CAinfo_storage[specindex] = 's'; // sparse matrix (sparse)
      else CAinfo_storage[specindex] = 'd'; // low-rank matrix (dense)
      if(blktype[blk-1] == 'd') CAinfo_type[specindex] = 'd'; // this should never occur in this input format
      else if(nnz[specindex] > 0) CAinfo_type[specindex] = 's';
      else CAinfo_type[specindex] = 'l';

    }

  CAinfo_entptr[(m+1)*numblk] = ct;
  CAinfo_rowcolptr[(m+1)*numblk] = ct;

  // Allocate remaining space plus some temp space for use in Step 3

  MYCALLOC(CAent, double, ct);
  MYCALLOC(CArow, size_t, ct);
  MYCALLOC(CAcol, size_t, ct);


  ////////////////////////////////
  // Step 3
  ////////////////////////////////

  // Open file
  datafile = fopen(datafilename, "r");

  // Redo first lines in file
  ret = fscanf(datafile, "%d", &m);
  ret = fscanf(datafile, "%d", &numblk);
  for(blk = 0; blk < numblk; blk++) {
    ret = fscanf(datafile, "%d", &(blksz[blk])); 
    if(blksz[blk] < 0) blksz[blk] *= -1;
  }
  for(i = 0; i < m; i++) ret = fscanf(datafile, "%lf", &(b[i]));
  ret = fscanf(datafile, "%lf", &entry);

  // Read entries in file size_to data structures
  // This needs work.
  for(i = 0; i <= m; i++) {
    for(k = 0; k < numblk; k++) {
      ret = fscanf(datafile, "%d %d %c %d", &num, &blk, &type, &sz);
      ct = CAinfo_entptr[DATABLOCKIND(i,blk,numblk)];
      if(type == 's') {
        for(j = 1; j <= sz; j++) {
          ret = fscanf(datafile, "%d %d %lf", &pos1, &pos2, &entry);
          if(fabs(entry) > DBL_EPSILON) {
            if(pos1 == 0 || pos2 == 0) { printf("Error (readdata_sdplr): Encountered '0' row or column index.\n"); exit(0); }
            if(pos1 > pos2) { CArow[ct] = pos1; CAcol[ct] = pos2; }
            else { CArow[ct] = pos2; CAcol[ct] = pos1; }
            CAent[ct] = entry;
            ct++;
          }
        }
      }
      else if(type == 'l') {
        for(j = 1; j <= sz*(blksz[blk-1] + 1); j++) {
          ret = fscanf(datafile, "%lf", &entry);
          if(j <= sz) { pos1 = j; pos2 = j; }
          else {
            pos1 = (j - sz)%blksz[blk-1];
            pos2 = (j - sz)/blksz[blk-1] + 1;
            if(pos1 == 0) { pos1 = blksz[blk-1]; pos2--; }
          }
          CArow[ct] = pos1; CAcol[ct] = pos2;
          CAent[ct] = entry;
          ct++;
        }
      }
    }
  }

  // Finished with this pass through file
  fclose(datafile);

  /////////////////////////////////
  // End Step 3
  /////////////////////////////////

  /////////////////////////////////
  // End Step 2
  /////////////////////////////////

  /////////////////////////////////
  // End Step 1
  /////////////////////////////////

  MYFREE(nnz);

  /////////////////////////////////
  // End Step 0
  /////////////////////////////////

  // Return back to calling function

  *passed_m                = m;
  *passed_numblk           = numblk;
  *passed_blksz            = blksz;
  *passed_blktype          = blktype;
  *passed_b                = b;
  *passed_CAent            = CAent;
  *passed_CArow            = CArow;
  *passed_CAcol            = CAcol;
  *passed_CAinfo_entptr    = CAinfo_entptr;
  *passed_CAinfo_rowcolptr = CAinfo_rowcolptr;
  *passed_CAinfo_type      = CAinfo_type;
  *passed_CAinfo_storage   = CAinfo_storage;

  return 0;
}


size_t readdata_raw(char* datafilename, size_t* passed_m, size_t* passed_numblk, ptrdiff_t** passed_blksz,
                char** passed_blktype, double** passed_b, double** passed_CAent,
                size_t** passed_CArow, size_t** passed_CAcol, size_t** passed_CAinfo_entptr,
                size_t** passed_CAinfo_rowcolptr, char** passed_CAinfo_type,
                char** passed_CAinfo_storage)
{
  size_t h, i, k;
  size_t ret;
  FILE *fid;

  size_t    m, numblk;
  ptrdiff_t *blksz;
  char   *blktype;
  double *b, *CAent;
  size_t    *CArow, *CAcol, *CAinfo_entptr, *CAinfo_rowcolptr;
  char   *CAinfo_type, *CAinfo_storage;

  fid = fopen(datafilename, "r");
  if(fid == NULL) {
    printf("Can't get file %s\n", datafilename);
    exit(0);
  }

  ret = fscanf(fid, "%d\n", &m);

  ret = fscanf(fid, "%d\n", &numblk);

  MYCALLOC(blksz, ptrdiff_t, numblk);
  MYCALLOC(blktype, char, numblk);
  MYCALLOC(b, double, m);

  for(k = 1; k <= numblk; k++)
    ret = fscanf(fid, "%d %c\n", &(blksz[k-1]), &(blktype[k-1]));

  for(i = 1; i <= m; i++)
    ret = fscanf(fid, "%lf\n", &(b[i-1]));

  MYCALLOC(CAinfo_entptr, size_t, (m+1)*numblk + 1);
  MYCALLOC(CAinfo_rowcolptr, size_t, (m+1)*numblk + 1);
  MYCALLOC(CAinfo_type, char, (m+1)*numblk);
  MYCALLOC(CAinfo_storage, char, (m+1)*numblk);

  for(h = 0; h < (m+1)*numblk; h++)
    ret = fscanf(fid, "%d %d %c %c\n", &(CAinfo_entptr[h]), &(CAinfo_rowcolptr[h]), &(CAinfo_type[h]), &(CAinfo_storage[h]));
  
  ret = fscanf(fid, "%d %d\n", &(CAinfo_rowcolptr[(m+1)*numblk]), &(CAinfo_entptr[(m+1)*numblk]));

  MYCALLOC(CArow, size_t, CAinfo_rowcolptr[(m+1)*numblk]);
  MYCALLOC(CAcol, size_t, CAinfo_rowcolptr[(m+1)*numblk]);
  MYCALLOC(CAent, double, CAinfo_entptr[(m+1)*numblk]);

  for(h = 0; h < CAinfo_rowcolptr[(m+1)*numblk]; h++)
    ret = fscanf(fid, "%d %d\n", &(CArow[h]), &(CAcol[h]));

  for(h = 0; h < CAinfo_entptr[(m+1)*numblk]; h++)
    ret = fscanf(fid, "%lf\n", &(CAent[h]));

  // Return back to calling function

  *passed_m                = m;
  *passed_numblk           = numblk;
  *passed_blksz            = blksz;
  *passed_blktype          = blktype;
  *passed_b                = b;
  *passed_CAent            = CAent;
  *passed_CArow            = CArow;
  *passed_CAcol            = CAcol;
  *passed_CAinfo_entptr    = CAinfo_entptr;
  *passed_CAinfo_rowcolptr = CAinfo_rowcolptr;
  *passed_CAinfo_type      = CAinfo_type;
  *passed_CAinfo_storage   = CAinfo_storage;

  return 0;

}


size_t writedata_raw(char* datafilename, size_t m, size_t numblk, ptrdiff_t* blksz,
                  char* blktype, double* b, double* CAent,
                  size_t* CArow, size_t* CAcol, size_t* CAinfo_entptr,
                  size_t* CAinfo_rowcolptr, char* CAinfo_type,
                  char* CAinfo_storage)
{
  size_t h, i, k;
  FILE *fid;

  fid = fopen(datafilename, "w");
  if(fid == NULL) {
    printf("Warning (writedata_raw): Could not open file for writing.\n");
    return 0;
  }

  fprintf(fid, "%d\n", m);

  fprintf(fid, "%d\n", numblk);

  for(k = 1; k <= numblk; k++)
    fprintf(fid, "%d %c\n", blksz[k-1], blktype[k-1]);

  for(i = 1; i <= m; i++)
    fprintf(fid, "%.16e\n", b[i-1]);

  for(h = 0; h < (m+1)*numblk; h++)
    fprintf(fid, "%d %d %c %c\n", CAinfo_entptr[h], CAinfo_rowcolptr[h], CAinfo_type[h], CAinfo_storage[h]);
  
  fprintf(fid, "%d %d\n", CAinfo_rowcolptr[(m+1)*numblk], CAinfo_entptr[(m+1)*numblk]);

  for(h = 0; h < CAinfo_rowcolptr[(m+1)*numblk]; h++)
    fprintf(fid, "%d %d\n", CArow[h], CAcol[h]);

  for(h = 0; h < CAinfo_entptr[(m+1)*numblk]; h++)
    fprintf(fid, "%.16e\n", CAent[h]);

  fclose(fid);

  return 0;
}

size_t writedata_sdplr(char* datafilename, size_t m, size_t numblk, ptrdiff_t* blksz,
                    char* blktype, double* b, double* CAent,
                    size_t* CArow, size_t* CAcol, size_t* CAinfo_entptr,
                    char* CAinfo_type)
{
  size_t h, i, j, k;
  size_t r, c, nnz, sz, info, lwork, rank=0, maxsz;
  char jobz, uplo;
  double *MM, *w, *work, maxeval=0.0;
  FILE *fid;

//   if(numblk != 1) {
//     printf("Error (writedata_sdplr): Format currently only supports one block.\n");
//     return 0;
//   }

  // Anything having to do with eigen assumes only one block
  
  maxsz = -1;
  for(k = 0; k < numblk; k++)
    if(blksz[k] > maxsz) maxsz = blksz[k];

  jobz = 'V';
  uplo = 'L';
  lwork = 3*maxsz-1;
  MYCALLOC(MM, double, maxsz*maxsz);
  MYCALLOC(w, double, maxsz);
  MYCALLOC(work, double, lwork);

  fid = fopen(datafilename, "w");
  if(fid == NULL) {
    printf("Warning (writedata_sdplr): Could not open file for writing.\n");
    return 0;
  }

  fprintf(fid, "%d\n", m);

  fprintf(fid, "%d\n", numblk);

  for(k = 1; k <= numblk; k++) {
    if(blktype[k-1] == 's') fprintf(fid, "%d\n", blksz[k-1]);
    else if(blktype[k-1] == 'd') fprintf(fid, "%d\n", -blksz[k-1]);

  }

  for(i = 1; i <= m; i++)
    fprintf(fid, "%.16e  ", b[i-1]);
  fprintf(fid, "\n");

  fprintf(fid, "-1.0\n");

  for(h = 0; h <= m; h++)
  {
    for(k = 1; k <= numblk; k++)
    {
      if(CAinfo_type[ DATABLOCKIND(h,k,numblk) ] == 's')
      {
        // Experimental: check for low rank if density is high
        sz  = blksz[k-1];
        nnz = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ];
        if((double)2*nnz/(sz*(sz+1)) - 0.75 > DBL_EPSILON) {
          for(j = 0; j < sz*sz; j++) MM[j] = 0.0;
          for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++) {
            r = CArow[j] - 1;
            c = CAcol[j] - 1;
            MM[c*sz + r] = CAent[j];
            MM[r*sz + c] = CAent[j];
          }
          dsyev_(&jobz, &uplo, &sz, MM, &sz, w, work, &lwork, &info);
          if(info == 0) {
            maxeval = -1.0e10;
            for(j = 0; j < sz; j++) if(fabs(w[j]) - maxeval > DBL_EPSILON) maxeval = fabs(w[j]);
            rank = 0;
            for(j = 0; j < sz; j++) if(fabs(w[j])/maxeval > DBL_EPSILON) rank++;
            printf("(h,k) = (%d,%d) : rank %d\n", h, k, rank);
          }
//           else printf("h = %d : eval computation bad\n", h);

        }

        if(rank <= sz/10 && (double)2*nnz/(sz*(sz+1)) - 0.75 > DBL_EPSILON) {
          fprintf(fid, "%d %d l %d\n", h, k, rank);
          for(j = 0; j < sz; j++) if(fabs(w[j]/maxeval) > DBL_EPSILON)
            fprintf(fid, "%.15e\n", w[j]);
          for(j = 0; j < sz; j++) if(fabs(w[j]/maxeval) > DBL_EPSILON)
            for(i = 0; i < sz; i++)
              fprintf(fid, "%.15e\n", MM[j*sz + i]);

        }
        else {
        
          // Right now, this matrix is assumed sparse
          j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ];
          fprintf(fid, "%d %d s %d\n", h, k, j);
          for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++)
            fprintf(fid, "%d %d %.16e\n", CArow[j], CAcol[j], CAent[j]);

        }

      }
      else if(CAinfo_type[ DATABLOCKIND(h,k,numblk) ] == 'l')
      {
        // Right now, this matrix is assumed dense
        j = (CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ])/(blksz[k-1]+1);
        fprintf(fid, "%d %d l %d\n", h, k, j);
        for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++)
          fprintf(fid, "%.16e\n", CAent[j]);
      }
      else if(CAinfo_type[ DATABLOCKIND(h,k,numblk) ] == 'd')
      {
        // Right now this matrix is assumed sparse
        j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ];
        fprintf(fid, "%d %d s %d\n", h, k, j);
        for(j = CAinfo_entptr[ DATABLOCKIND(h,k,numblk) ]; j <= CAinfo_entptr[ DATABLOCKIND(h,k,numblk) + 1 ] - 1; j++)
          fprintf(fid, "%d %d %.16e\n", CArow[j], CAcol[j], CAent[j]);
      }
      else {
        printf("Error (writedata_sdplr): Encountered data matrix not of type 's' or 'l' or 'd'.\n");
        fclose(fid);
        return 0;
      }
    }
  }

  fclose(fid);

  MYFREE(MM);
  MYFREE(w);
  MYFREE(work);

  return 0;

}

size_t writedata_sdpa(char* datafilename, size_t m, size_t numblk, ptrdiff_t* blksz,
                   char* blktype, double* b, double* CAent,
                   size_t* CArow, size_t* CAcol, size_t* CAinfo_entptr,
                   char* CAinfo_type)
{
  size_t h, i, j, k;
  FILE *fid;
  int tmp;

  printf("writedata_sdpa: Warning! The sense of the optimization may be wrong.\n");

  fid = fopen(datafilename, "w");
  if(fid == NULL) {
    printf("Warning (writedata_sdpa): Could not open file for writing.\n");
    return 0;
  }

  fprintf(fid, "%d\n", m);

  fprintf(fid, "%d\n", numblk);

  for(k = 1; k <= numblk; k++) {
    if(blktype[k-1] == 's') fprintf(fid, "%d ", blksz[k-1]);
    else if(blktype[k-1] == 'd') fprintf(fid, "%d ", -blksz[k-1]);
  }
  fprintf(fid, "\n");

  for(i = 1; i <= m; i++)
/*     fprintf(fid, "%.16e  ", b[i-1]); */
    fprintf(fid, "%.0f  ", b[i-1]);
  fprintf(fid, "\n");

  for(h = 0; h <= m; h++)
  {
    for(k = 1; k <= numblk; k++)
    { 
      tmp = DATABLOCKIND(h,k,numblk);
      if(CAinfo_type[tmp] == 's')
      {
        // Right now, this matrix is assumed sparse
        for(j = CAinfo_entptr[tmp]; (int)j <= (int)CAinfo_entptr[tmp+1] - 1; j++) {
          if(h == 0) fprintf(fid, "%d %d %d %d %.16e\n", h, k, CArow[j], CAcol[j], -CAent[j]);
          else       fprintf(fid, "%d %d %d %d %.16e\n", h, k, CArow[j], CAcol[j],  CAent[j]);
        }
      }
      else if(CAinfo_type[tmp] == 'l')
      {
        printf("error: Low-rank matrices not supported in SDPA format.\n");
        exit(0);
      }
      else if(CAinfo_type[tmp] == 'd')
      {
        // Right now this matrix is assumed sparse
        for(j = CAinfo_entptr[tmp]; (int)j <= (int)CAinfo_entptr[tmp+1] - 1; j++)
          if(h == 0) fprintf(fid, "%d %d %d %d %.16e\n", h, k, CArow[j], CAcol[j], -CAent[j]);
          else fprintf(fid, "%d %d %d %d %.16e\n", h, k, CArow[j], CAcol[j], CAent[j]);
      }
      else {
        printf("Error (writedata_sdplr): Encountered data matrix not of type 's' or 'l' or 'd'.\n");
        fclose(fid);
        return 0;
      }
    }
  }

  fclose(fid);

  return 0;

}
