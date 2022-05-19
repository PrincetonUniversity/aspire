// Compute the correlations between all rays of two projections, without 
// taking shifts into account.
//
// This function has no test code, as this kernel was replaced by it 
// successors, such as kernel6.cu and on.
//
// Yoel Shkolnisky, May 2022.

// Add two complex numbers.
// Given two complex numbers a and b, reuturb a+b.     
__device__ __forceinline__ float2 ComplexAdd(float2 a, float2 b)
{
    float2 c;

    c.x = a.x + b.x;
    c.y = a.y + b.y;

    return c;
}

// Multiply two complex numbers.
// Given two complex numbers a and b, reuturb a*b.     
__device__ __forceinline__ float2 ComplexMul(float2 a, float2 b)
{
    float2 c;

    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;

    return c;
}

// Multiply two complex numbers with conjugation.
// Given two complex numbers a and b, reuturb conj(a)*b.
__device__ __forceinline__ float2 ComplexMulConj(float2 a, float2 b)
{
    float2 c;

    c.x = a.x * b.x + a.y * b.y;
    c.y = a.x * b.y - a.y * b.x;

    return c;
}

// Multiply two matrices
// Given A of size nxp, B of size pxm, and a target matrix C of size nxm, 
// compute C=A*B.
__device__ __forceinline__ void ComplexMatMul(const float2* A, const float2* B, 
        float2* C, const int n, const int m, const int p){
// A of size nxp, B of size pxm, C of size nxm.
    float2 ab;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < p; ++k){
                ab = ComplexMul(A[i+n*k],B[k+p*j]);
                C[i+n*j] = ComplexAdd(C[i+n*j],ab);
            }
}

// Multiply two matrices while conjugating the first matrix.
// Given A of size pxn, B of size pxm, and a target matrix C of size nxm, 
// compute C=A'*B. (prime is conjugate transpose).
__device__ __forceinline__ void ComplexMatMulConj(const float2* A, const float2* B, 
        float2* C, const int n, const int m, const int p){
// A of size pxn, B of size pxm, C of size nxm.
    float2 ab;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < p; ++k){
                ab = ComplexMulConj(A[k+p*i],B[k+p*j]);
                C[i+n*j] = ComplexAdd(C[i+n*j],ab);
            }
}

// Multiply two matrices while transposing (without conjugation) the first matrix.
// Given A of size pxn, B of size pxm, and a target matrix C of size nxm, 
// compute C=A'*B. (prime is conjugate transpose).
__device__ __forceinline__ void ComplexMatMulTrans(const float2* A, const float2* B, 
        float2* C, const int n, const int m, const int p){
// A of size pxn, B of size pxm, C of size nxm.
    float2 ab;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            for (int k = 0; k < p; ++k){
                ab = ComplexMul(A[k+p*i],B[k+p*j]);
                C[i+n*j] = ComplexAdd(C[i+n*j],ab);
            }
}


// Take the real part of a matrix.
// A is complex with n elements. reA is float with n elements.
// Note that if A and reA are matrices, then n should be the total number
// of elements in the matrices.        
__device__ __forceinline__ void real(const float2* A, float* reA, const int n){
    for (int i=0; i<n; i++)
        reA[i] = A[i].x;
}

// Pointwise multiplication of a matrix by vector. The vector is expanded
// to all columns of the matrix.
// A is of size nxm. V is of size n. AV is of size nxm and is equal to 
// AV=bsxfun(@times,A,V);        
__device__ __forceinline__ void MatVecPointwise(const float2* A, const float2* V, 
        float2* AV, const int n, const int m){
    for (int j=0; j<m; j++){
        for (int i=0; i<n; i++){
            AV[i+j*n] = ComplexMul(A[i+j*n],V[i]);
        }
    }
}

// Find the maximum of a float array A with n elements.
// Returns the maximal value and the index of the maximum.
// If the same maximal value occurs several times, the first occurrence is 
// returned.         
// Note that the output index variable is of type single, so that all types
// in the calling function are of type single (to avoid mixing int32).        
__device__ __forceinline__ void max(const float* A, const int n, float* val, float* idx ){
    *val=A[0];
    *idx=0.0f;

    for (int i=1; i<n; i++){
        if (A[i]>*val){
            *val=A[i];
            *idx=(float) i;
        }
    }
}

__global__ void commonline(const float2* P1, const float2* P2, 
        const int n_r, const int n_theta, float *C1, float* C2, 
        float* cij, float* cji, float* corrstack){
    
    float ab;
    float maxcorr = -2.0; // Maximal correlation observed thus far.
    float cl12,cl21;        // The common line index at which the maxiaml 
            // correlation was obtained. These are floats so it will be 
            // easier to return to the calling host.
    float corr; // Correlation of the currently inpsected pair of lines.
    
// C1=2*real(P1_shifted'*P2);
// Essentially implemeting  
//        2.*real(ComplexMatMulConj(P1, P2, C2tmp, n_theta,n_theta,n_r));
// without allocating an addition array        
         
    for (int i = 0; i < n_theta; ++i){
        for (int j = 0; j < n_theta; ++j){
            corr = 0.0;
            for (int k = 0; k < n_r; ++k){
                // We only compute the real part of the product
                ab = P1[k+n_r*i].x * P2[k+n_r*j].x + P1[k+n_r*i].y * P2[k+n_r*j].y;
                corr += ab;
            }
            corr *= 2.0;
            C1[i+n_theta*j] = corr;

            if (corr>maxcorr){ // Current pair of lines is a better candidate for the common line
                    maxcorr = corr;
                    cl12 = (float) i;
                    cl21 = (float) j;
            }
        }
    }

//  C2=2*real(P1_shifted_flipped'*P2);
// Essentially implemeting  
//        2.*real(ComplexMatMulTrans(P1, P2, C2tmp, n_theta,n_theta,n_r));
// without allocating an addition array        
    for (int i = 0; i < n_theta; ++i){
        for (int j = 0; j < n_theta; ++j){
            corr = 0.0;
            for (int k = 0; k < n_r; ++k){
                // We only compute the real part of the product
                ab = P1[k+n_r*i].x * P2[k+n_r*j].x - P1[k+n_r*i].y * P2[k+n_r*j].y;
                corr += ab;
            }
            corr *= 2.0;
            C2[i+n_theta*j] = corr;

            if (corr>maxcorr){ // Current pair of lines is a better candidate for the common line
                    maxcorr = corr;
                    cl12 = (float) i;
                    cl21 = (float) j+n_theta;
            }
        }
    }        
    // Return results
    *cij = cl12 + 1.0;  // Host indexing is 1-based.
    *cji = cl21 + 1.0;
    *corrstack  = maxcorr;
}
