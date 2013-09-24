#include "myinclude.h"
#ifdef __MEX
#include "mex.h"
#endif

#define NUMPARAMS 14

// Default params (still hard coded in generate_params)
#define INPUTTYPE     1
#define RHO_F         1.0e-5
#define RHO_C         1.0e-1
#define SIGMAFAC      2.0
#define RANKREDUCE    0
#define TIMELIM       3600
#define PRINTLEVEL    1
#define DTHRESH_DIM   10
#define DTHRESH_DENS  0.75
#define NUMBFGSVECS   4
#define RANKREDTOL    2.2204460492503131e-16
#define GAPTOL        1.0e-3
#define CHECKBD       -1
#define TYPEBD        1

// size_t generate_params_info(size_t);
// size_t getparams_maxlinelength(FILE*);
// size_t getparams_getline(FILE*, char*, size_t);
// size_t getparams_tolower(char*, size_t);
// size_t getparams_isnumber(char*);

// Macros just for this file
// #define MYCALLOC(VAR,TYPE,SIZE) VAR = (TYPE*)calloc(SIZE, sizeof(TYPE))
// #define MYFREE(VAR) free(VAR)

size_t getparams(
char*  paramfile,
size_t*    inputtype,
double* rho_f,
double* rho_c,
double* sigmafac,
size_t*    rankreduce,
size_t*    timelim,
size_t*    printlevel,
size_t*    dthresh_dim,
double* dthresh_dens,
size_t*    numbfgsvecs,
double* rankredtol,
double* gaptol,
ptrdiff_t*    checkbd,
size_t*    typebd)
{
  size_t i, buffsz, ret, numparams=NUMPARAMS;
  size_t assigned[NUMPARAMS];
  double values[NUMPARAMS];
  char *buff, *match;
  FILE *fid;

  char paramstr[NUMPARAMS][50] = {
    "input type",
    "feasibility tolerance",
    "centrality tolerance",
    "penalty factor increase",
    "rank reduction",
    "time limit",
    "print level",
    "dimension threshold for dense matrices",
    "density threshold for dense matrices",
    "number of lbfgs vectors",
    "rank tolerance",
    "duality gap tolerance",
    "check dual bound",
    "dual bound type"
  };

  // Assign default params just to be sure.
  *inputtype    = INPUTTYPE;         values[0]  = (double)INPUTTYPE; 
  *rho_f        = RHO_F;             values[1]  = (double)RHO_F;     
  *rho_c        = RHO_C;             values[2]  = (double)RHO_C;     
  *sigmafac     = SIGMAFAC;          values[3]  = (double)SIGMAFAC;  
  *rankreduce   = RANKREDUCE;        values[4]  = (double)RANKREDUCE;
  *timelim      = TIMELIM;           values[5]  = (double)TIMELIM;   
  *printlevel   = PRINTLEVEL;        values[6] =  (double)PRINTLEVEL;
  *dthresh_dim  = DTHRESH_DIM;       values[7]  = (double)DTHRESH_DIM;
  *dthresh_dens = DTHRESH_DENS;      values[8]  = (double)DTHRESH_DENS;
  *numbfgsvecs  = NUMBFGSVECS;       values[9]  = (double)NUMBFGSVECS;
  *rankredtol   = RANKREDTOL;        values[10] = (double)RANKREDTOL;
  *gaptol       = GAPTOL;            values[11] = (double)GAPTOL; 
  *checkbd      = CHECKBD;           values[12] = (double)CHECKBD;
  *typebd       = TYPEBD;            values[13] = (double)TYPEBD;


  // Simple case (param file = NULL)
  if(paramfile == NULL) return 1;

  // Keep track of what has been assigned by file.
  for(i = 0; i < numparams; i++) assigned[i] = 0;

  // Open file for reading
  fid = fopen(paramfile, "r");

  if(fid == NULL) {
    printf("Warning (getparams): File %s not found. Using default parameters.\n", paramfile);
    return 0;
  }
  else if(fid != NULL) {

    buffsz = getparams_maxlinelength(fid) + 10;
    fclose(fid);

    fid = fopen(paramfile, "r");

    MYCALLOC(buff, char, buffsz);

    do {
      ret = getparams_getline(fid, buff, buffsz);
      getparams_tolower(buff, buffsz);

      for(i = 0; i < numparams; i++) {
        match = strstr(buff, paramstr[i]);
        if(match != NULL) {
          if(assigned[i] == 0) {
            match = strchr(buff, ':');
            if(match == NULL) {
              printf("Error (getparams): Parameter file has wrong format.\n");
              return -1;
            }
            match++;
            if(getparams_isnumber(match) != 1) {
              printf("Error (getparams): Parameter file has wrong format.\n");
              return -1;
            }
            values[i] = atof(match);
            assigned[i] = 1;
          }
          else if(assigned[i] == 1)
            printf("Warning (getparams): Attempt to assign parameter '%s' twice.\n", paramstr[i]);
        }
      }

    } while(ret != 0);

    MYFREE(buff);

    fclose(fid);

    // If some values have not been assigned, print warning
    for(i = 0; i < numparams; i++)
      if(assigned[i] == 0)
        printf("Warning (getparams): Some parameters not assigned. Using default values.\n");

    // assign values read from file
    *inputtype    = (size_t)values[0];
    *rho_f        =      values[1];
    *rho_c        =      values[2];
    *sigmafac     =      values[3];
    *rankreduce   = (size_t)values[4];
    *timelim      = (size_t)values[5];
    *printlevel   = (size_t)values[6];
    *dthresh_dim  = (size_t)values[7];
    *dthresh_dens =      values[8];
    *numbfgsvecs  = (size_t)values[9];
    *rankredtol   =      values[10];
    *gaptol       =      values[11];
    *checkbd      = (ptrdiff_t)values[12];
    *typebd       = (size_t)values[13];

    // Perform some simple checks

    if(*inputtype != 1 && *inputtype != 2 && *inputtype != 1000) {
      printf("Error (params): Parameter '%s' must be 1 or 2.\n", paramstr[0]);
      return -1;
    }

    if(*rho_f <= 0) {
      printf("Error (params): Parameter '%s' must be positive.\n", paramstr[1]);
      return -1;
    }

    if(*rho_c <= 0) {
      printf("Error (params): Parameter '%s' must be positive.\n", paramstr[2]);
      return -1;
    }

    if(*sigmafac <= 1.0) {
      printf("Error (params): Parameter '%s' must be greater than 1.0.\n", paramstr[4]);
      return -1;
    }

    if(*rankreduce != 0 && *rankreduce != 1) {
      printf("Error (params): Parameter '%s' must be 0 or 1.\n", paramstr[5]);
      return -1;
    }

    if(*timelim <= 0) {
      printf("Parameter '%s' must be positive.\n", paramstr[5]);
      return -1;
    }

    if(*printlevel != 0 && *printlevel != 1) {
      printf("Error (params): Parameter '%s' must be 0 or 1.\n", paramstr[6]);
      return -1;
    }

    if(*dthresh_dim < 0) {
      printf("Parameter '%s' must be nonnegative.\n", paramstr[7]);
      return -1;
    }

    if(-(*dthresh_dens) > DBL_EPSILON || *dthresh_dens - 1.0 > DBL_EPSILON) {
      printf("Parameter '%s' must be in [0,1].\n", paramstr[8]);
      return -1;
    }

    if(*numbfgsvecs < 0) {
      printf("Error (params): Parameter '%s' must be a non-negative size_teger.\n", paramstr[9]);
      return -1;
    }

    if(*rankredtol <= 0.0) {
      printf("Error (params): Parameter '%s' must be positive.\n", paramstr[10]);
      return -1;
    }

    if(*gaptol <= 0) {
      printf("Error (params): Parameter '%s' must be positive.\n", paramstr[11]);
      return -1;
    }

    // Hopefully temporary
    if(*checkbd != -1) {
      printf("Error (params): At this time, parameter '%s' must be -1.\n", paramstr[12]);
      return -1;
    }

    if(*checkbd != -1 && *checkbd != 0 && *checkbd != 1) {
      printf("Error (params): Parameter '%s' must be -1, 0, or 1.\n", paramstr[12]);
      return -1;
    }

    if(*typebd != 1) {
      printf("Error (params): Currently, parameter '%s' must equal 1.\n", paramstr[13]);
      return -1;
    }

    return 1;

  }

  return -1;

}

size_t getparams_maxlinelength(FILE *datafile)
{
  size_t maxlen=0, k, c;
  
  do {

    k = 0;

    do {
      c = getc(datafile); 
      k++;
    } while(c != '\n' && c != EOF);

    if(k > maxlen) maxlen=k;

  } while(c != EOF);
  
  return(maxlen);
}

size_t getparams_getline(FILE *datafile, char *buffer, size_t bufsiz)
{
  size_t k;
  char c;
  
  k = 0;

  do {
    if(k >= bufsiz) {
      printf("Error (getparams_getline): Line too long!  Adjust bufsiz.\n");
      return(-1);
    }
    c = getc(datafile);
    buffer[k] = c;
    k++;
  } while(c != '\n' && c != EOF);

  if(c != '\n') {
    buffer[k] = '\n';
    buffer[k+1] = 0;
  }
  else buffer[k] = 0;

  if(c == EOF) return 0;
  else return 1;
}

size_t getparams_tolower(char* buff, size_t buffsz)
{
  size_t i;

  for(i = 0; i < buffsz; i++)
    buff[i] = tolower(buff[i]);

  return 1;
}

size_t getparams_isnumber(char* str)
{
  size_t i, length;

  length = strlen(str);

  for(i = 0; i < length; i++)
    if(isdigit(str[i]) == 0 && str[i] != '.' &&
       str[i] != '-' && str[i] != 'e' && isspace(str[i]) == 0 &&
       str[i] != '\n' && str[i] != '\0' && str[i] != EOF && str[i] != '+') {
      return 0;
    }

  return 1;
}

size_t generate_params(void)
{

  /*
  size_t     inputtype;
  double  rho_f;
  double  rho_c;
  double  sigmafac;
  size_t     rankreduce;
  size_t     timelim;
  size_t     printlevel;
  size_t     dthresh_dim;
  double  dthresh_dens;
  size_t     numbfgsvecs;
  double  rankredtol;
  double  gaptol;
  ptrdiff_t     checkbd;
  size_t     typebd;
  */

  size_t i;
  char in[NUMPARAMS][110], *ret, filename[] = "sdplr.params", usrfilename[100];
  FILE *fid;
  char paramstr[NUMPARAMS][100] = {
    "Input type (1=SDPA, 2=SDPLR)",
    "Feasibility tolerance",
    "Centrality tolerance",
    "Penalty factor increase",
    "Rank reduction (1=yes, 0=no)",
    "Time limit (in seconds)",
    "Print level (0=nothing, 1=all)",
    "Dimension threshold for dense matrices",
    "Density threshold for dense matrices (in [0,1])",
    "Number of LBFGS vectors",
    "Rank tolerance (e.g., machine precision)",
    "Duality gap tolerance",
    "Check dual bound (-1, 0, 1)",
    "Dual bound type (1=+/-1, 2=0/1)"
  };
  char def[NUMPARAMS][50] = {
    "1",
    "1.0e-5",
    "1.0e-1",
    "2.0",
    "0",
    "3600",
    "1",
    "10",
    "0.75",
    "4",
    "2.2204460492503131e-16",
    "1.0e-3",
    "-1",
    "1"
  };

  printf("\nSDPLR %s  --  Automatic Paramater File Generation\n\n", VERSION);

  do {
    printf("\n");
    printf("Parameter file name [%s]: ", filename); fflush(stdout);
    ret = fgets(usrfilename, 100, stdin);
    if(ret == NULL) {
      printf("Error\n");
      exit(0);
    }
    usrfilename[strlen(usrfilename)-1] = '\0';
    if(strlen(usrfilename) == 0) strcpy(usrfilename, filename);
    fid = fopen(usrfilename, "w");
  } while(fid == NULL);

  printf("\n\nPress 'i' for information at any time.\n");
  printf("Press 'i' for information at any time.\n");
  printf("Press 'i' for information at any time.\n\n"); fflush(stdout);

  for(i = 0; i < NUMPARAMS; i++) {
    do {
      printf("\n");
      printf("%s [%s]: ", paramstr[i], def[i]); fflush(stdout);
      ret = fgets(in[i], 100, stdin);
      if(ret == NULL) {
        printf("Error\n");
        exit(0);
      }
      in[i][strlen(in[i])-1] = '\0';
      if(strlen(in[i]) == 0) strcpy(in[i], def[i]);
      if(in[i][0] == 'i' || in[i][0] == 'I') generate_params_info(i);
    } while(getparams_isnumber(in[i]) == 0);
  }

  fprintf(fid, "SDPLR %s paramter file (automatically generated)\n\n", VERSION);
  fprintf(fid, "--> Basic parameters <--\n\n");
  for(i = 0; i < 10; i++)
    fprintf(fid, "%s : %s\n", paramstr[i], in[i]);
  fprintf(fid, "\n--> In-development parameters <--\n\n");
  for(i = 10; i < NUMPARAMS; i++)
    fprintf(fid, "%s : %s\n", paramstr[i], in[i]);

  fclose(fid);

  printf("\n");

  return 0;

}

size_t generate_params_info(size_t num)
{
  switch(num)
  {
  case 0:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies which format your datafile is in, either\n");
      printf("SDPA or SDPLR format. SDPA is the most common format.\n\n");
      printf("Information on the SDPA format can be found at\n");
      printf("http://www.nmt.edu/~sdplib/FORMAT. Information on the SDPLR format\n");
      printf("can be found in the SDPLR User's Guide.\n");
      printf("\n");
      break;
    }
  case 1:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies how accurately you would like to satisfy\n");
      printf("the primal constrasize_ts. SDPLR is a primal infeasible method that\n");
      printf("works towards feasibility as the algorithm progresses. Primal\n");
      printf("infeasibility is measured in relative terms as\n");
      printf("\n");
      printf("            || b - A(RR^T) || / (1 + |b_max|)\n");
      printf("\n");
      printf("where R is the decision variable, A represents the constrasize_ts, b\n");
      printf("is the right-hand side, and b_max is the largest entry of b in\n");
      printf("absolute value.\n");
      printf("\n");
      printf("Smaller values will cause SDPLR to do more work. On most problems,\n");
      printf("reasonable values are in the range 1.0e-5 to 1.0e-8.\n");
      printf("\n");
      break;
    }
  case 2:
    {
      printf("\n"); printf("\n");
      printf("This parameter specifies how accurately you would like to solve\n");
      printf("each augmented Lagrangian subproblem. SDPLR uses the augmented\n");
      printf("Lagrangian approach for nonlinear problems. Optimality for a\n");
      printf("subproblem is measured in relative terms as\n");
      printf("\n");
      printf("            || gradient || = || 2SR ||_F / (1 + |C_max|)\n");
      printf("\n");
      printf("where R is the decision variable, S is the dual variable estimate,\n");
      printf("C is the objective cost matrix, and C_max is the largest entry of\n");
      printf("C in absolute value. Based on the above formula, this parameter\n");
      printf("can be size_terpreted as enforcing complementary slackness between\n");
      printf("the primal and dual problems.\n");
      printf("\n");
      printf("Smaller values will cause SDPLR to do more work. On most problems,\n");
      printf("reasonable values are in the range 1.0e-1 and 1.0e-3.\n");
      printf("\n");
      break;
    }
  case 3:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the factor by which the penalty parameter\n");
      printf("is increased every major iteration. SDPLR uses the penalty\n");
      printf("parameter to enforce feasibility in its augmented Lagrangian\n");
      printf("approach. This parameter should be greater than 1.0.\n");
      printf("\n");
      printf("Smaller values are considered more conservative; higher values are\n");
      printf("more aggressive. Reasonable values are between 2.0 and 10.0.\n");
      printf("\n");
      break;
    }
  case 4:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies whether or not to perform the rank\n");
      printf("reduction procedure, which is a dynamic way to reduce the\n");
      printf("dimensionality of the problem, thereby speeding up the algorithm.\n");
      printf("On many problems, rank reduction is very effective, but on some\n");
      printf("problems it can actually increase the overall time required.\n");
      printf("\n");
      break;
    }
  case 5:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the overall time limit for the algorithm\n");
      printf("(in seconds). The algorithm is terminated immediately after the\n");
      printf("completion of the first minor iteration past the time limit.\n");
      printf("\n");
      break;
    }
  case 6:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the amount of information printed to the\n");
      printf("screen by the algorithm: 0 prints nothing, and 1 prints\n");
      printf("everything.\n");
      printf("\n");
      break;
    }
  case 7:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the threshold N such that all matrices\n");
      printf("with both dimensions <= N are stored as dense matrices. The idea\n");
      printf("is that, even if a small matrix is sparse, it is often quicker to\n");
      printf("treat it as dense.\n");
      printf("\n");
      printf("In the current implementation of SDPLR, this parameter applies\n");
      printf("only to storage of the dual variable S, and not the data matrices\n");
      printf("C, A_i, which are always stored as sparse. The primal variable R\n");
      printf("is always stored as dense.\n");
      printf("\n");
      break;
    }
  case 8:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the threshold d such that all matrices\n");
      printf("with density >= d are stored as dense matrices. Here, d is in\n");
      printf("[0,1].\n");
      printf("\n");
      printf("In the current implementation of SDPLR, this parameter applies\n");
      printf("only to storage of the dual variable S, and not the data matrices\n");
      printf("C, A_i, which are always stored as sparse. The primal variable R\n");
      printf("is always stored as dense.\n");
      printf("\n");
      break;
    }
  case 9:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies how many limited memory BFGS (LBFGS)\n");
      printf("pairs to store. A higher number produces better LBFGS directions\n");
      printf("but requires more computation. Reasonable values are between 3 and\n");
      printf("10.\n");
      printf("\n");
      break;
    }
  case 10:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* Only applicable if parameter 'rank reduction' is set to\n");
      printf("yes=1.\n");
      printf("\n");
      printf("This parameter controls the determination of numerical rank of the\n");
      printf("primal variable R in the rank reduction procedure. Let the\n");
      printf("parameter be called TOL. In essence, if a certain submatrix M (k\n");
      printf("columns) of R (r columns) satisfies\n");
      printf("\n");
      printf("           || M ||_F <= TOL * || R ||_F,\n");
      printf("\n");
      printf("then the numerical rank of R is r-k. A good choice for TOL is the\n");
      printf("machine precision. If TOL is set too high, then the rank reduction\n");
      printf("procedure may become unreliable.\n");
      printf("\n");
      break;
    }
  case 11:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* This parameter refers to a feature of SDPLR that is in\n");
      printf("development but currently inactive! It has to do with calculating\n");
      printf("dual bounds with SDPLR. Contact samuel-burer@uiowa.edu for more\n");
      printf("information.\n");
      printf("\n");
      printf("This parameter specifies the target duality gap. SDPLR terminates\n");
      printf("once both the target feasibility and gap are met.\n");
      printf("\n");
      break;
    }
  case 12:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* This parameter refers to a feature of SDPLR that is in\n");
      printf("development but currently inactive! It has to do with calculating\n");
      printf("dual bounds with SDPLR. Contact samuel-burer@uiowa.edu for more\n");
      printf("information.\n");
      printf("\n");
      printf("This parameter specifies when to calculate a dual bound. The value\n");
      printf("-1 means never, 0 means only at the end of the algorithm, 1 means\n");
      printf("after every major iteration.\n");
      printf("\n");
      break;
    }
  case 13:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* This parameter refers to a feature of SDPLR that is in\n");
      printf("development but currently inactive! It has to do with calculating\n");
      printf("dual bounds with SDPLR. Contact samuel-burer@uiowa.edu for more\n");
      printf("information.\n");
      printf("\n");
      printf("This parameter specifies what type of dual bound to calculate. The\n");
      printf("value 1 refers to +/-1 combinatorial bounds, and the value 2\n");
      printf("refers to 0/1 combinatorial bounds.\n");
      printf("\n");
      break;
    }
  default:
    {
      printf("default\n");
      break;
    }
  }

  fflush(stdout);

  return 0;

}



