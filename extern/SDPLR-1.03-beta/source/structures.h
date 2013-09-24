//! \brief Core iteration information for implementation of L-BFGS.
//!
//! SDPLR uses the limited memory BFGS (L-BFGS) method as its core
//! unconstrained minimization procedure. This method is based on
//! storing several types of vectors and scalars for each iteration of
//! the algorithm. These vectors and scalars are more or less linear
//! combinations of posize_ts and previous posize_ts. This structure is the
//! basic collection of information for one iteration. In SDPLR, an
//! array of structures of this type will be masize_tained and updated.
//! This size of this array is \ref NUMBFGSVECS, which is the so-called
//! "memory" of the L-BFGS procedure and is a user-specified size_teger
//! parameter.

// #ifndef __MEX
// #define size_t int
// #endif

typedef struct {

  //! \brief Difference between iterates
  //!
  //! This stores the difference between successive iterates, i.e.,
  //! \f$ x^{k+1} - x^k \f$. (Here, \f$ x \f$ is generic.)

  double *s; 

  //! \brief Difference between gradients
  //!
  //! This stores the difference between successive gradients, i.e.,
  //! \f$ \nabla f(x^{k+1}) - \nabla f(x^k) \f$. (Here, \f$ x \f$ and
  //! \f$ f \f$ are generic.)
 
  double *y;

  //! \brief Reciprocal of inner product between \ref s and \ref y
  //!
  //! This stores \f$ (y^T s)^{-1} \f$. It is significant only as
  //! a temporary storage places.

  double rho;

  //! \brief Temporary storage.
  //!
  //! Significant only as a temporary storage place related to \ref s,
  //! \ref y, and \ref rho.
  
  double a;

} lbfgsvec;

typedef struct {
  double*  d;
  double*  ent;
  size_t      nrow;
  size_t      ncol;
} lowrankmat;

typedef struct {
  size_t*    row;
  size_t*    col;
  size_t     nnz;
  double* ent;
  size_t*    XS_in;
} sparsesymmmat;

typedef struct {
  size_t*    ind;
  size_t     nnz;
  double* ent;
  size_t*    XS_in;
} diagmat;

typedef struct {
  lowrankmat*    lr;
  sparsesymmmat* sp;
  diagmat*       diag;
  char           type;
  double*        multval;
  char*          label;
} datamat;

typedef struct {

  // user options
  double     rho_f;
  double     rho_c;
  double     sigmafac;
  size_t        rankreduce;
  size_t        timelim;
  size_t        printlevel;
  size_t        dthresh_dim;
  double     dthresh_dens;
  size_t        numbfgsvecs;
  double     rankredtol;
  double     gaptol;
  size_t        checkbd;
  size_t        typebd;

  // very basic data (both SDPA and lowrank formats)
  size_t        m;
  size_t        numblk;
  size_t*       blksz;
  char*      blktype;
  datamat*** A;
  double*    b;
  datamat**  C;

  // auxiliary data, applies to lowrank format only
  double     evaladj;

  //
  // algorithm data, doesn't change once it is setup
  //
  // global dimension
  size_t        n;
  size_t        nr;
  // number of limited memory vecs to store
  size_t        M;
  // posize_ters to different types of matrices
  size_t**      lrind;
  size_t*       nlrind;
  size_t**      usind;
  size_t*       nusind;
  // rank information
  size_t*       rank;
  size_t*       maxrank;
  // norm of C for size_ternal stopping criterion
  double     normC;


  //
  // algorithm data, changes per iteration
  //
  // basic structures
  double*      lambda;
  double       sigma;
  double*      vio;
  double*      G;

  // timing data structures
  double totaltime;
  clock_t start;
  clock_t finish;

  // What are these?
  double   *S;
  double   *X;
  double   *y;
  size_t      *XS_blkptr;
  char     *XS_blksto;
  size_t     **XS_colptr;
  size_t     **XS_rowind;
  size_t      *AA_rowptr;
  size_t      *AA_colind;
  double   *AA_colval_one;
  double   *AA_colval_two;
  size_t      *lr_mat;
  size_t      *lr_blk;
  size_t       lr_num;

} problemdata;



