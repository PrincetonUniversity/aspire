// Play with writing CU files that contain not only a single __global__ 
// function, but also several __device__ functions that can be called.
// Yoel Shkolnisky, May 2022.
        
        
// Add two complex numbers.
// Given two complex numbers a and b, reuturb a+b.     
__device__ inline float2 ComplexAdd(float2 a, float2 b)
{
    float2 c;

    c.x = a.x + b.x;
    c.y = a.y + b.y;

    return c;
}

// Multiply two complex numbers.
// Given two complex numbers a and b, reuturb a*b.     
__device__ inline float2 ComplexMul(float2 a, float2 b)
{
    float2 c;

    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;

    return c;
}

// Multiply two complex numbers with conjugation.
// Given two complex numbers a and b, reuturb conj(a)*b.
__device__ inline float2 ComplexMulConj(float2 a, float2 b)
{
    float2 c;

    c.x = a.x * b.x + a.y * b.y;
    c.y = a.x * b.y - a.y * b.x;

    return c;
}

// Multiply two matrices
// Given A of size nxp, B of size pxm, and a target matrix C of size nxm, 
// compute C=A*B.
__device__ inline void ComplexMatMul(const float2* A, const float2* B, 
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
__device__ inline void ComplexMatMulConj(const float2* A, const float2* B, 
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

// Take the real part of a matrix.
// A is complex with n elements. reA is float with n elements.
// Note that if A and reA are matrices, then n should be the total number
// of elements in the matrices.        
__device__ inline void real(const float2* A, float* reA, const int n){
    for (int i=0; i<n; i++)
        reA[i] = A[i].x;
}

// Pointwise multiplication of a matrix by vector. The vector is expanded
// to all columns of the matrix.
// A is of size nxm. V is of size n. AV is of size nxm and is equal to 
// AV=bsxfun(@times,A,V);        
__device__ inline void MatVecPointwise(const float2* A, const float2* V, 
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
__device__ inline void max(const float* A, const int n, float* val, float* idx ){
    *val=A[0];
    *idx=0.0f;

    for (int i=1; i<n; i++){
        if (A[i]>*val){
            *val=A[i];
            *idx=(float) i;
        }
    }
}

/*
__global__ void clmatrix(const float2* A, const float2* B, float2* C, 
        const int n, const int m, const int p){
            ComplexMatMul(A, B, C, n, m, p);
}
__global__ void myreal(const float2* A, float* reA, int n){
    real(A,reA,n);
}

__global__ void bsxfun(const float2* A, const float2* V, float2* AV, const int n, const int m){
    MatVecPointwise(A, V, AV, n, m);
}

*/
        
__global__ void callmax(const float* A, const int n, float* val, float* idx){
    max(A,n,val,idx);
}
        