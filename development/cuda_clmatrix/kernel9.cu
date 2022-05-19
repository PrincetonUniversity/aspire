// Like kernel8.cu but with further optimizations.
// This kernel is called by cryo_clmatrix_gpu_kernel_3.m, which is tested 
// by test_cryo_clmatrix_gpu_kernel_3.m
//
// Yoel Shkolnisky, May 2022.
        

#define PI_F 3.141592654f

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

// Exponential of a complex number.
__device__ __forceinline__ float2 ComplexExp(float2 z)

{
    float2 res;
    float t = expf (z.x);
    sincosf (z.y, &res.y, &res.x);
    res.x *= t;
    res.y *= t;
    return res;
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

__global__ void commonline(const float2* P, const int n_projs,
        const int n_r, const int n_theta, const float max_shift, 
        const float shift_step, float* clstack, float* corrstack){
    
    int k1 =  blockIdx.x * blockDim.x + threadIdx.x;
    if (k1 >= n_projs) return;

    int k2 = blockIdx.y * blockDim.y + threadIdx.y;
    if (k2 >= n_projs) return;

    if (k1>=k2) return;

    const float2 *P1=P+k1*n_r*n_theta;
    const float2 *P2=P+k2*n_r*n_theta;
            
            
    float ab1,ab2;
    float maxcorr = -2.0; // Maximal correlation observed thus far.
    float cl12,cl21;        // The common line index at which the maxiaml 
            // correlation was obtained. These are floats so it will be 
            // easier to return to the calling host.
    float corr1,corr2; // Correlation of the currently inpsected pair of lines.
             
    float n_shifts = ceilf(2.0*max_shift/shift_step+1.0); // Number of shifts to try.
            
    for (int shiftidx=0; shiftidx<n_shifts; shiftidx++){
        float shift=-max_shift+shiftidx*shift_step;

        //shift_phases=exp(-2*pi*sqrt(-1).*rk2.*shift./(2*rmax+1)); 
        float2 c0 = ComplexExp(make_float2(0,-2.0*PI_F*(-n_r)*shift/(2*n_r+1)));
        float2 c1 = ComplexExp(make_float2(0,-2.0*PI_F*shift/(2*n_r+1)));


        for (int i = 0; i < n_theta; ++i){
            for (int j = 0; j < n_theta; ++j){
                        
                corr1 = 0.0;
                corr2 = 0.0;
                float2 phi = c0;

                for (int k = 0; k < n_r; ++k){

                    float2 P1_shifted;
                    float2 P1ik = P1[k+n_r*i];
                    float2 P2jk = P2[k+n_r*j];
                    // C1=2*real(P1_shifted'*P2);
                    // Essentially implemeting  
                    //        2.*real(ComplexMatMulConj(P1, P2, C2tmp, n_theta,n_theta,n_r));
                    // without allocating an addition array        
                            
                    // We only compute the real part of the product
                    P1_shifted = ComplexMul(phi,P1ik);
                    ab1 = P1_shifted.x * P2jk.x + P1_shifted.y * P2jk.y;
                    corr1 += ab1;

                    //  C2=2*real(P1_shifted_flipped'*P2);
                    // Essentially implemeting  
                    //        2.*real(ComplexMatMulTrans(P1, P2, C2tmp, n_theta,n_theta,n_r));
                    // without allocating an addition array        

                    P1_shifted = ComplexMulConj(phi,P1ik);
                    ab2 = P1_shifted.x * P2jk.x - P1_shifted.y * P2jk.y;
                    corr2 += ab2;

                    phi = ComplexMul(phi,c1); // Generate next shift phase
                }
                corr1 *= 2.0;
                corr2 *= 2.0;


                float s,a,b;
                s = copysignf(1.0,corr1-maxcorr);
                a = (1.0+s)/2.0;
                b = (1.0-s)/2.0;
    
                maxcorr = a*corr1 + b*maxcorr;
                cl12 = a*i + b*cl12;
                cl21 = a*j + b*cl21;

                s = copysignf(1.0,corr2-maxcorr);
                a = (1.0+s)/2.0;
                b = (1.0-s)/2.0;
    
                maxcorr = a*corr2 + b*maxcorr;
                cl12 = a*i + b*cl12;
                cl21 = a*(j+n_theta) + b*cl21;

            }
        }

    }
    // Return results
    clstack[k1+k2*n_projs] = cl12 + 1.0;  // Host indexing is 1-based.
    clstack[k2+k1*n_projs] = cl21 + 1.0;
    corrstack[k1+k2*n_projs]  = maxcorr;  // Writing to global memory is very slow. So return only what is necessary.
}
    