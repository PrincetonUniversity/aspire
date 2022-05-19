// Kernel implementation of clmatrix_cheat
//
// Yoel Shkolnisky, May 2022

#define PI 3.1415926f
        
#define rmulAB(A,B,C)\
        C[0] = A[6]*B[2] + A[3]*B[1] + A[0]*B[0];\
        C[1] = A[7]*B[2] + A[4]*B[1] + A[1]*B[0];\
        C[2] = A[8]*B[2] + A[5]*B[1] + A[2]*B[0];\
        C[3] = A[6]*B[5] + A[3]*B[4] + A[0]*B[3];\
        C[4] = A[7]*B[5] + A[4]*B[4] + A[1]*B[3];\
        C[5] = A[8]*B[5] + A[5]*B[4] + A[2]*B[3];\
        C[6] = A[6]*B[8] + A[3]*B[7] + A[0]*B[6];\
        C[7] = A[7]*B[8] + A[4]*B[7] + A[1]*B[6];\
        C[8] = A[8]*B[8] + A[5]*B[7] + A[2]*B[6];
        
#define mulABT(A,B,C)\
        C[0] = A[6]*B[6] + A[3]*B[3] + A[0]*B[0];\
        C[1] = A[7]*B[6] + A[4]*B[3] + A[1]*B[0];\
        C[2] = A[8]*B[6] + A[5]*B[3] + A[2]*B[0];\
        C[3] = A[6]*B[7] + A[3]*B[4] + A[0]*B[1];\
        C[4] = A[7]*B[7] + A[4]*B[4] + A[1]*B[1];\
        C[5] = A[8]*B[7] + A[5]*B[4] + A[2]*B[1];\
        C[6] = A[6]*B[8] + A[3]*B[5] + A[0]*B[2];\
        C[7] = A[7]*B[8] + A[4]*B[5] + A[1]*B[2];\
        C[8] = A[8]*B[8] + A[5]*B[5] + A[2]*B[2];

__global__ void clmatrix_cheat(const float* rots, int n_rots, int L, float* clmatrix){
                
                
    // k1 is the index of rots1 to process.    
    int k1 =  blockIdx.x * blockDim.x + threadIdx.x;
    if (k1 >= n_rots-1) return;

    // k2 is the index of rots2 to process.
    int k2 = blockIdx.y * blockDim.y + threadIdx.y;
    if (k2 >= n_rots) return;

    // We only process the indices idx=0:n_rots-2, k2=k1+1:N-1
    if (k2 <= k1) return;

    // Get points to the processed rotations.
    const float *Ri = rots+k1*9;
    const float *Rj = rots+k2*9;
      
    float Ut[9];

    // Ut=Rj*Ri.';
    mulABT(Rj,Ri,Ut)
            

    // alphaij=atan2(Ut(3,1),-Ut(3,2));
    // alphaji=atan2(-Ut(1,3),Ut(2,3));
    float alphaij = atan2(Ut[2],-Ut[5]);
    float alphaji = atan2(-Ut[6],Ut[7]);

    // alphaij=alphaij+PI; % Shift from [-pi,pi] to [0,2*pi].
    // alphaji=alphaji+PI;
    alphaij=alphaij+PI;
    alphaji=alphaji+PI;

    // l_ij=alphaij/(2*PI)*L;
    // l_ji=alphaji/(2*PI)*L;

    float l_ij=alphaij/(2*PI)*L;
    float l_ji=alphaji/(2*PI)*L;

    // l_ij=mod(round(l_ij),L);
    // l_ji=mod(round(l_ji),L);

    l_ij=fmodf(round(l_ij),L);
    l_ji=fmodf(round(l_ji),L);
                        
    // Save output
    clmatrix[k2*n_rots + k1] = l_ij+1;
    clmatrix[k1*n_rots + k2] = l_ji+1;
}