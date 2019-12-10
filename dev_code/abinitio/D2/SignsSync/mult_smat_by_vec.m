
% mult_smat_by_vec applies the signs sync matrix to a vector using: 
% Input: rows_arr= an array of size 3NxN where rows_arr(3*i-2:3*i,j) is a
%        row of rotation Ri calculated using iterative scheme with image j
%        pairs_map= array of size NxN where
%        pairs_map(i,j)==uppertri_ijtoind(i,j)

function [v_out]=mult_smat_by_vec(v,sign_mat,pairs_map,N)

v_out=zeros(size(v));

for i=1:N
    for j=i+1:N    
       ij=uppertri_ijtoind(i,j,N);
       v_out(ij)=sign_mat(ij,:)*v(pairs_map(ij,:));
    end
end

