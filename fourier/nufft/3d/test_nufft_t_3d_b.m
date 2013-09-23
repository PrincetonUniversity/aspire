function test_nufft_t_3d_b
% 
% Compare the results of the non-equally spaced FFT to the three
% dimensional FFT.
%
% Yoel Shkolnisky, January 2009.
%
% Revisions:
% Filename changed from test_nufft_t3_2 to test_nufft_t_3d_b. The call to
% nufft_t3_v2 was changed to a call to nufft_t_3d. (Y.S. December 22, 2009)
%

n=32;
v=rand(n,n,n);
range=-fix(n/2):n-fix(n/2)-1;
[J1,J2,J3]=ndgrid(range,range,range);
J=-2*pi*[J1(:) J2(:) J3(:)]/n;

f1=nufft_t_3d(v,J,'single');
f1=reshape(f1,n,n,n);
f2=cfftn(v);

disp(norm(f1(:)-f2(:))/norm(f2(:)))

