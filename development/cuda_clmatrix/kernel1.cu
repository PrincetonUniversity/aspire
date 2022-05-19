// Make sure we know how to access comples numbers.
// Get a complex 2D array and return two real arrays of the same size 
// containing the real and imaginary parts.        
//
// Yoel Shkolnisky, May 2022
        
__global__ void corr(const float2* comp_in, int n_r, int n_theta, float* real_out, float* imag_out ) {
    // comp_in is a a complex array of size n_rxn_teta. 
    // The function returns two arrays of the same size containing its real 
    // and imaginary parts.
    
    int ridx =  blockIdx.x * blockDim.x + threadIdx.x;
    if (ridx >= n_r) return;

    int cidx = blockIdx.y * blockDim.y + threadIdx.y;
    if (cidx >= n_theta) return;

    real_out[cidx*n_r+ridx]= (float) comp_in[cidx*n_r+ridx].x;
    imag_out[cidx*n_r+ridx]= (float) comp_in[cidx*n_r+ridx].y;
}