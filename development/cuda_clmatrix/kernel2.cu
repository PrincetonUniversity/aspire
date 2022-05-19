// Try loops in CUDA kernel.
// Compute all inner products between the comlumns of the two complex arrays.
//
// Yoel Shkolnisky, May 2022
        
__global__ void corr(const float2* pf1, const float2* pf2, int n_r, int n_theta, float2* corrs) {
    // pf1 and pf2 are two complex arrays of size n_r X n_thets. 
    // corrs is a complex array of size n_theta x n_theta whose (i,j) element
    // is the inner product of pf1(:,i) and pf2(:,j).
            
    // idx1 is the index of the column in pf1 to process.    
    int idx1 =  blockIdx.x * blockDim.x + threadIdx.x;
    if (idx1 >= n_theta) return;

    // idx2 is the index of the column in pf2 to process.
    int idx2 = blockIdx.y * blockDim.y + threadIdx.y;
    if (idx2 >= n_theta) return;

    int i1 = idx1*n_r; // Base address in pf1
    int i2 = idx2*n_r; // Base address in pf2
            
    float2 c = make_float2(0.0,0.0);
    for (int i=0; i<n_r; i++){
        c.x += pf1[i1].x * pf2[i2].x + pf1[i1].y * pf2[i2].y;
        c.y += pf1[i1].y * pf2[i2].x - pf1[i1].x * pf2[i2].y;
        i1++; i2++; 
    }
 
    corrs[idx2*n_theta + idx1] = c;
}