# This Makefile.inc works on Ubuntu 9.04 (April 2009)

################################
# LAPACK directory and library #
################################

# If a relative path is given next, it should be relative from the
# source/ subdirectory.

LAPACK_LIB_DIR=
LAPACK_LIB= -llapack

##############################
# BLAS directory and library #
##############################

# If a relative path is given next, it should be relative from the
# source/ subdirectory.

BLAS_LIB_DIR=
BLAS_LIB= -lblas

####################
# Main definitions #
####################

CC       = gcc
CFLAGS   = -O3
LIB_DIRS = -L$(LAPACK_LIB_DIR) -L$(BLAS_LIB_DIR)
LIBS     = -lgsl $(LAPACK_LIB) $(BLAS_LIB) -lgfortran -lm
