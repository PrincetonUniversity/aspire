
%Changed to zeros(...,'single') 06/01/2019
function [cMat_eff]=fill_signSync_matrix_c(cMat5Dmp,sync_signs2,img,c,N)
cMat_eff=zeros(3,3,N,N);
for r=1:N %row of matrix to fill
    cMat_eff(:,:,r,:)=cMat5Dmp(:,:,r,sync_signs2(r,:,img,c),c);
end