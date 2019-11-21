function [ CR ] = rotateCellFast( cell1, angles, Wp, Wm )

% Fast rotation of using 'Pinchon Vector' for coefficients from spherical 
% expansion of a volume
%
% INPUT: 
%   cell1: cell1{ll+1,1} represents the l-th order coefficient where m 
%       stands for column index, and radial index n stands for row index
%   angles: 3*N vector representing Euler ZYZ angles of desired rotations
%   Wp: wigner D matrix for Euler ZYZ angle [-pi/2, pi/2, pi/2], 
%   Wm: wigner D matrix for Euler ZYZ angle [-pi/2, -pi/2, pi/2],
%
% OUTPUT: 
%   CR: spherical expansion after rotation, of size 

maxL = size(cell1,1)-1;
N = size(angles,2);

CR = cell(maxL+1,1);

for ll = 0:maxL
    atemp = reshape(transpose(cell1{ll+1,1}),2*ll+1,1,[]);
    atemp = bsxfun(@times, exp(-1i*(-ll:ll)'*angles(1,:)), atemp);
    atemp = reshape(atemp, 2*ll+1,[]);
    atemp = Wp{1,ll+1}*atemp;
    atemp = reshape(atemp, 2*ll+1, N, []);
    atemp = bsxfun(@times, exp(-1i*(-ll:ll)'*angles(2,:)), atemp);
    atemp = reshape(atemp, 2*ll+1, []);
    atemp = Wm{1,ll+1}*atemp;
    atemp = reshape(atemp, 2*ll+1, N, []);
    atemp = bsxfun(@times, exp(-1i*(-ll:ll)'*angles(3,:)), atemp);
    CR{ll+1,1} = permute(atemp(ll+1:end,:,:), [3,1,2]);
end

end

