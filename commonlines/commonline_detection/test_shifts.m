function [ est_shifts, err] = test_shifts( shift_equations, ref_shifts)
% Test the accuracy of the shift_equations/estimated shifts.
% 
% Input: 
%       shift_equations  System of equations for the 2D shifts of the
%       projections. This is a sparse system with 2*n_proj+1 columns. The
%       first 2*n_proj columns correspond to the unknown shifts dx,dy of
%       each projection. The last column is the right-hand-side of the
%       system, which is the relative shift of each pair of common-lines.
%
%       ref_shifts: true shifts.
%
% Output:
%       est_shifts: estimated shifts.
%       err: l2 norm errors of the estimated shifts.
%
% Lanhui Wang, July 2, 2013
%
K = size(ref_shifts,1); % number of images

% Shift estimation using LS, and a l2 regularization on sizes of shifts.
% Three additional equations on first three shift parameters fix the 3
% degreees of freedom of the linear system.
est_shifts=[shift_equations(:,1:end-1);sparse(1:3,1:3,ones(1,3),3,2*K);eye(2*K)]...
    \[shift_equations(:,end);ref_shifts(1,1);ref_shifts(1,2);ref_shifts(2,1);zeros(2*K,1)];
est_shifts = reshape(est_shifts, 2, K)';
err = sqrt(sum((est_shifts - ref_shifts).^2,2));
end

