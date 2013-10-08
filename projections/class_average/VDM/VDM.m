function [ class_VDM, class_VDM_refl, angle ] = VDM(class, corr, rot, class_refl, k, flag, n_nbor) 
%Function for using VDM to search for nearest neighbors more robust to
%noise
% Input:
%       class: 
%           Pxl matrix, where l>k. This gives a list of l-nearest
%           neighbors for P data points.
%       corr: 
%           Pxl matrix, normalized cross correlation between each data
%           point and its l-nearest neighbors.
%       rot: 
%           Pxl matrix, in-plane rotational alignment in degrees between
%           each data point and its l-nearest neighbors.
%       class_refl: 
%           indicates whether reflection is included. class_refl==1: no
%           reflection. class_refl==2: reflection.
%       k: 
%           number of nearest neighbors to choose to make the graph.
%       flag:
%           indicates generating k-nn graph using union rule (flag==0) or
%           k-nn graph using linear programming.
%       n_nbor: 
%           number of nearest neighbors user needs.
%Output:
%       class_VDM: 
%           Pxn_nbor matrix. It gives n_nbor nearest neighbors for each
%           data point.
%       class_VDM_refl:
%           Pxn_nbor matrix. It indicates whether there is a reflection.
%           class_VDM_refl==1: no reflection. class_VDM_refl==2:
%           reflection.
%       angle:
%           Pxn_nbor matrix. The rotational alignment in degrees between
%           each point and its nearest neighbors.
%
%Zhizhen Zhao Aug 2013


P=size(class, 1); 

[ X, Ah, rows, cols ] = sort_list_weights_wrefl(class(:, 1:k), sqrt(2-2*real(corr(:, 1:k))), rot(:, 1:k), class_refl(:, 1:k) );
if flag==0
    W = ones(length(X), 1); %Generating knn graph using union rule
else
    W =script_find_graph_weights_v3(X, [rows, cols], 2*P, 5); %Generating k-nn graph using linear programming. It makes sure that the weight in each row is the same.
end;
W2=W.*Ah;
H2=sparse(rows(W>0.001), cols(W>0.001), W2(W>0.001), 2*P, 2*P);
H2=H2+H2';

[r_VDM_LP, ~, VV_LP]=VDM_LP(H2, 24);

if P<=10^4
    corr_VDM=r_VDM_LP(1:P, :)*r_VDM_LP';
    corr_VDM=real(corr_VDM-sparse(1:P, 1:P, ones(P, 1), P, 2*P));
    [~, class_VDM]=sort(corr_VDM, 2, 'descend');
    class_VDM=class_VDM(1:P, 1:n_nbor);
else
    P_max=5000;
    for i=1:ceil(P/P_max)
        corr_VDM=r_VDM_LP((i-1)*P_max+1:min(i*P_max, P), :)*r_VDM_LP';
        corr_VDM=real(corr_VDM);
        [~, tmp]=sort(corr_VDM, 2, 'descend');
        class_VDM((i-1)*P_max+1:min(i*P_max, P), :)=tmp(:, 2:n_nbor+1);
    end;
end;
        
class_VDM_refl=ceil(class_VDM/P);
class_VDM(class_VDM>P)=class_VDM(class_VDM>P)-P;
list=[repmat([1:P]', n_nbor, 1), class_VDM(:)+(class_VDM_refl(:)-1)*P];
[ angle ] = VDM_angle_v2(VV_LP(:, 1:10), list);

angle=reshape(angle, P, n_nbor);

end

