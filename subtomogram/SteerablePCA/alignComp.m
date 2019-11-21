function [eul1 ] = alignComp( vec1, vec2, siz_n )
% Align two eigenvectors representing the principal component in expansion
% coefficient.

% Input
% vec1: coefficient to align
% vec2: coefficient for reference volume

% Output
% rotated_vec1: rotated coefficient of vec1 according to vec2

c1 = coeff2mat(cell2coeff(vec2cell(vec1, siz_n)));
c2 = coeff2mat(cell2coeff(vec2cell(vec2, siz_n)));
L = size(siz_n,1);
[ func_2d, ~] = SOFT_yuan( c1, c2, L, 20*L);
[~,ind] = max(real(func_2d));
eul1 = ind2euler(ind, L);


end

