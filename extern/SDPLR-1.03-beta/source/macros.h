#ifndef __WIN32
#define _isnan isnan
#endif

#ifdef __MEX
#define printf mexPrintf
#endif

// #ifndef __MEX
// #define size_t int
// #endif

#define MYCALLOC(VAR,TYPE,SIZE) VAR = (TYPE*)calloc(SIZE, sizeof(TYPE))
#define MYFREE(VAR) free(VAR);

#define mymaxint(A, B) ((A) > (B) ? (A) : (B))
// #define mymin(A, B) ((A) < (B) ? (A) : (B))
#define mymax(A, B) ((A) - (B) > DBL_EPSILON ? (A) : (B))
#define mymin(A, B) ((B) - (A) > DBL_EPSILON ? (A) : (B))

#define SMATIND(I,J,DIM) (J-1)*DIM + I
#define LMATIND(I,J,DIM) DIM*(DIM+1)/2 -(DIM-J+1)*(DIM-J+2)/2 + I - J + 1

#define EASYDAXPY(DIM,SCAL,X,Y)   mydaxpy(DIM,SCAL,X+1,1,Y+1,1)
#define EASYDCOPY(DIM,X,Y)        mydcopy(DIM,X+1,1,Y+1,1)
#define EASYDDOT(DIM,X,Y)         myddot(DIM,X+1,1,Y+1,1)
#define EASYDNRM2(DIM,X)          mydnrm2(DIM,X+1,1)
#define EASYDSCAL(DIM,SCAL,X)     mydscal(DIM,SCAL,X+1,1)
#define ZEROVEC(X,DIM)            mydscal(DIM,0.0,X+1,1)

// Some parameters for testing purposes
#define RANDOM 0 // default 1 (also in main.c (???))
#define SCALE_OBJ 0

#define SDPBLK  's'
#define DIAGBLK 'd'
#define SPARSE  's'
#define DENSE   'd'

#define DATABLOCKIND(DATA,BLOCK,NUMBLOCK) ((DATA+1)-1)*NUMBLOCK + BLOCK - 1

#ifdef __WIN32
#ifndef __MINGW32__
#define dsyr_ dsyr
#define dsyrk_ dsyrk
#define dsyr2k_ dsyr2k
#define dgemm_ dgemm
#define dsymm_ dsymm
#define idamax_ idamax
#define dsyev_ dsyev
#define dgeqp3_ dgeqp3
#define daxpy_ daxpy
#define dcopy_ dcopy
#define ddot_ ddot
#define dnrm2_ dnrm2
#define dscal_ dscal
#endif
#endif


#ifdef __MEX
#undef exit
#define exit(er) mexErrMsgTxt("SDPLR Error\n");
#endif
