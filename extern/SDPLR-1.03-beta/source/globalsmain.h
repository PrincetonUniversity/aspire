// global variables

//! \brief Global variable that keeps track of SDP block. Only necessary
//! if \ref __ARPACK is defined.
//!
//! Arpack requires a function posize_ter to complete its MAT*vec
//! operation, but as far as I know, there's no good way to pass
//! parameters size_to that function (something I would definitely like to
//! do). So my work-around is to make global copies of those parameters
//! and access them directly inside Arpack. See \ref my_arpack and \ref
//! simple_Stimesvec_block.

// #ifndef __MEX
// #define size_t int
// #endif

size_t           global_blk;

//! \brief Global variable that is a simple copy of \ref problemdata
//! structure used elsewhere. Only necessary if \ref __ARPACK is
//! defined.
//!
//! See \ref global_blk for full explanation.

problemdata  *global_data;

//! \brief Global variable that is used as temporary space during data
//! operation \ref Aoper of \ref dataoper.c. Specifically used to help
//! with low-rank data matrices.
//!
//! The function \ref Aoper requires temporary space for operations such
//! as \f$ RR^T \f$ and \f$ RD^T \f$, where \f$ R \f$ is the variable
//! and \f$ D \f$ is the direction. This global variable provides that
//! space, and since it is global, it only has to be allocated once.
//! This is in hopes of improving performance as \ref Aoper is called
//! many times.

double       *global_UtB;

//! \brief Global variable that is used as temporary space during data
//! operation \ref Aoper of \ref dataoper.c. Specifically used to help
//! with low-rank data matrices.
//!
//! See \ref global_UtB for full explanation.

double       *global_VtB;

//! \brief Global variable that is used as temporary space during the
//! exact linesearch of \ref linesearch.c.
//!
//! During \ref linesearch, functions such as \ref Aoper require
//! temporary space for operations such as \f$ RR^T \f$ and \f$ RD^T \f$,
//! where \f$ R \f$ is the variable and \f$ D \f$ is the direction. This
//! global variable provides that space, and since it is global, it only
//! has to be allocated once. This is in hopes of improving performance
//! as \ref linesearch and \ref Aoper are called many times.
//! 
//! The size of this variable is related to the size of the XS
//! structures in \ref problemdata.

double *global_UVt;

//! \brief Global variable that is used as temporary space during the
//! exact linesearch of \ref linesearch.c.
//!
//! See \ref global_UVt for full explanation.

double *global_ARD;

//! \brief Global variable that is used as temporary space during the
//! exact linesearch of \ref linesearch.c.
//!
//! See \ref global_UVt for full explanation.

double *global_ADD;
