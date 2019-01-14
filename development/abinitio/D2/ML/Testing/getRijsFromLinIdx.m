
function [Rijs_est]=getRijsFromLinIdx(lookup_data,lin_idx)

Rijs_est=zeros(3,3,4,length(lin_idx));
l=size(lookup_data.cls_lookup,1)/4;
n_est_in_oct1=sum(lin_idx<=l);
if n_est_in_oct1>0
    Rijs_est11=getRijsFromLinIdxOct1(lookup_data,lin_idx(lin_idx<=l));
    Rijs_est(:,:,:,lin_idx<=l)=Rijs_est11;
end
if n_est_in_oct1<length(lin_idx)
    Rijs_est12=getRijsFromLinIdxOct12(lookup_data,lin_idx(lin_idx>l)-l);
    Rijs_est(:,:,:,lin_idx>l)=Rijs_est12;
end



