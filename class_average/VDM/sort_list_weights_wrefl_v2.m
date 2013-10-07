function [ X, Ah, rows, cols ] = sort_list_weights_wrefl_v2(class, corr, rot, refl )
%intersection rule
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
[~,ia, ib]=intersect(C, CC, 'rows');
rows=C(ia, 1);
cols=C(ia, 2);
rot_matrix=rot(:);
rot_matrix=[rot_matrix(ia)];
X=[X(ia)];

ind = find(rows<cols);
rows=rows(ind);
cols=cols(ind);
rot_matrix=rot_matrix(ind);
A=(rot_matrix)*2*pi/n_theta;
Ah=exp(sqrt(-1)*A);
X=X(ind);

end