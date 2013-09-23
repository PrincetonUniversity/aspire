/***********************************************************
 * Functions: printarr and printiarr                
 *    printarr  Print the matrix A of doubles of size mxn. 
 ***********************************************************/

#include <stdio.h>

void printarr(FILE* logfid, const double *A, int m, int n)
{
    int j,k;

    if (m*n>10000){
        fprintf(logfid,"Array too large. Not printint...\n");
        return;
    }
    
    for (j=0;j<m;j++){
        for (k=0;k<n;k++) {
            fprintf(logfid,"%+7.5E   ",A[k*m+j]);
        }
        fprintf(logfid,"\n");
    }
}


void printiarr(FILE* logfid, const int *A, int m, int n)
{
    int j,k;

    if (m*n>10000){
        fprintf(logfid,"Array too large. Not printint...\n");
        return;
    }

    for (j=0;j<m;j++){
        for (k=0;k<n;k++) {
            fprintf(logfid,"%3d   ",A[k*m+j]);
        }
        fprintf(logfid,"\n");
    }
}

