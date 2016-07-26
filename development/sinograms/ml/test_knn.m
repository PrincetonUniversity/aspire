rng(9823);
X=randn(100,1000); % 5000 signals in R^100
Y=randn(100,1000); % 5000 signals in R^100
k=10;

[N_ref,D_ref] = knn_opt(X, Y, k);
[N1,D1] = knn_opt_gpu_lowmem(X, Y, k,0.5e9);
[N2,D2]=  knn_opt_gpu_low_mem_2(X, Y, k,0.5e9);

assert((norm(D_ref(:)-D1(:))/norm(D_ref(:)))<1.0e-7);
assert(all(N_ref(:)-N1(:)==0))

assert((norm(D_ref(:)-D2(:))/norm(D_ref(:)))<1.0e-7);
assert(all(N_ref(:)-N2(:)==0));

% sort and mink break ties differently.
% Therefore, N may be slightly different from knn_opt and the other
% functions. 
% To compare N use
% for k=1:360; if ~isempty(setdiff(N(1:knn,k),N_ref(:,k))); disp(k); end;end
%
% Note that the differences in the D matrices are not zero but are very
% small (to single precision)