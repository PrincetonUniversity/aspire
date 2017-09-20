%TEST_ROT2XYZ Test the function test_rot2xyz.
%   The printed error should be of the order of machine precision.
%
%Yoel Shkolnisky, November 2013.
initstate;
maxerr=-1;
for k=1:10000
    R=rand_rots(1);
    [psi,theta,phi]=rot2xyz(R);
    R2=Rz(phi)*Ry(theta)*Rx(psi);
    err=norm(R2.'*R-eye(3),'fro');
    if err>maxerr
        maxerr=err;
    end
    if err>1
        aaa=1;
    end
end
fprintf('max error: %e\n',maxerr);
