function [ X, Ah, rows, cols ] = sort_list_weights_wrefl(class, corr, rot, refl )
%Function for computing graph structure.

n_theta=360;
P=size(class, 1);
n_nbor=size(class, 2);
class(refl==2)=class(refl==2)+P;
list=[class(:), repmat([1:P]', n_nbor, 1)];
class=class(:);
class(refl==1)=class(refl==1)+P;
class(refl==2)=class(refl==2)-P;
list_r = [ class(:), repmat([P+1:2*P]', n_nbor, 1)];
list=[list; list_r]; %including reflection
X=[corr(:); corr(:)]; %including reflection
rot=[rot(:); -rot(:)]; %including reflection
C=list;
CC=[C(:, 2), C(:, 1)];
[ia, ib] = union_row_idx(C, CC);
rows=[C(ia, 1); CC(ib, 1)];
cols=[C(ia, 2); CC(ib, 2)];
rot_matrix=rot(:);
rot_matrix=[rot_matrix(ia); -rot_matrix(ib)];
X=[X(ia); X(ib)];

ind = find(rows<cols);
rows=rows(ind);
cols=cols(ind);
rot_matrix=rot_matrix(ind);
A=(rot_matrix)*2*pi/n_theta;
Ah=exp(sqrt(-1)*A);
X=X(ind);

end

function [ia, ib] = union_row_idx(a, b)
    % Need to write this function ourselves since Octave (<4.2.0) has a bug
    % in its implementation.
    [~, idx] = unique([a; b], 'rows');

    a_rows = size(a, 1);

    ia = idx(idx <= a_rows);
    ib = idx(idx > a_rows) - a_rows;
end
