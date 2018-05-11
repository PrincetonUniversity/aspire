% Test speed and accuracy of optimizied voting functions.
%
% Test the optimized voting functions by comparing them to the old
% implementations.
%
% Yoel Shkolnisky, March 2015.

initstate;
TOL=1.0e-14;

Ks=[100, 200, 500, 1000];

for ntest=1:numel(Ks)
    K=Ks(ntest);
    rots_ref = rand_rots(K);
    L=360;     % Use a large number of lines per image, so we don't have discretization errors.
    cl=clmatrix_cheat(rots_ref,L);
    tic; 
    Sold=cryo_syncmatrix_vote_old(cl,L); 
    Told=toc;
    
    tic; 
    Snew=cryo_syncmatrix_vote(cl,L);
    Tnew=toc;
    
    speedup=Told/Tnew;
    err=norm(Sold-Snew)/norm(Sold);
    fprintf('Testing K=%d,\t Told=%d \t Tnew=%d \t speedup=%5.3f \t err=%7.5e\n',K,Told,Tnew,speedup,err);
    
end
