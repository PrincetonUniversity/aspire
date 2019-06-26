
function [oct1_ij_map,oct2_ij_map]=genSclsScoresIdxMap_eqClass_Dn(scls_lookup_data)

ntheta=scls_lookup_data.ntheta;
n=scls_lookup_data.n;
%eq_idx=scls_lookup_data.eq_idx;
%eq_idx2=scls_lookup_data.eq_idx2;
%eq_idx=single(eq_idx);
%eq_idx2=single(eq_idx2);
nrot=length(scls_lookup_data.scls_lookup1)/((2*n-1)*ntheta);
nrot2=length(scls_lookup_data.scls_lookup2)/((2*n-1)*ntheta);

%First the map of for i<j pairs from octatnt 1
%eq2eq_Rij_table=logical(triu(ones(nrot,nrot)-eq_idx*eq_idx',1));
eq2eq_Rij_table=scls_lookup_data.eq2eq_Rij_table11;
npairs=sum(eq2eq_Rij_table(:)); %pairs of candidates
oct1_ij_map=zeros(0.5*ntheta*ntheta,2,npairs);
i_idx=repmat(1:ntheta,0.5*ntheta,1);
j_idx=repmat((1:0.5*ntheta)',1,ntheta);
i_idx=i_idx(:);
j_idx=j_idx(:);
idx_vec=1:nrot;
idx=0;

for i=1:nrot
    
    unique_pairs_i=idx_vec(eq2eq_Rij_table(i,:));
    n2=sum(unique_pairs_i);
    if n2==0
        continue;
    end    
    i_idx_plus_offset=i_idx+(i-1)*ntheta;
    
    for j=unique_pairs_i   
        idx=idx+1;
        j_idx_plus_offset=j_idx+(j-1)*ntheta;
        oct1_ij_map(:,:,idx)=[i_idx_plus_offset,j_idx_plus_offset];
    end
end

%DEBUG 
% idx_vec=1:(0.5*ntheta^2*2*npairs);
% [j_theta,i_theta,ij,p]=ind2sub([0.5*ntheta,ntheta,2,npairs],idx_vec);


%Now construct the map for i in octant 1 with j in octant 2
%eq2eq_Rij_table=logical(ones(nrot,nrot2)-eq_idx*eq_idx2');
eq2eq_Rij_table=scls_lookup_data.eq2eq_Rij_table12;
npairs12=sum(eq2eq_Rij_table(:));
oct2_ij_map=zeros(0.5*ntheta*ntheta,2,npairs12);
idx_vec=1:nrot2;
idx=0;

for i=1:nrot
    
    unique_pairs_i=idx_vec(eq2eq_Rij_table(i,:));
    n2=sum(unique_pairs_i);
    if n2==0
        continue;
    end    
    i_idx_plus_offset=i_idx+(i-1)*ntheta;
    
    for j=unique_pairs_i   
        idx=idx+1;
        j_idx_plus_offset=j_idx+(j-1)*ntheta;
        oct2_ij_map(:,:,idx)=[i_idx_plus_offset,j_idx_plus_offset];
    end
end

tmp1=squeeze(oct1_ij_map(:,1,:));
tmp2=squeeze(oct1_ij_map(:,2,:));
oct1_ij_map=[tmp1(:),tmp2(:)];
tmp1=squeeze(oct2_ij_map(:,1,:));
tmp2=squeeze(oct2_ij_map(:,2,:));
oct2_ij_map=[tmp1(:),tmp2(:)];


