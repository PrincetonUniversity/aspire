

function [ij_map,ij_offset]=removeD2Ambiguities(sphere_grid)

%take care of first octant, consider all permutations of x,y,z no g_alpha
nproj=size(sphere_grid,2);
tol=1;%0.97;
ij_map=eye(nproj);
p=perms(3:-1:1);
p=p([4,5],:);
n_perms=size(p,1);

for k=1:n_perms
    gram_k=triu(sphere_grid'*sphere_grid(p(k,:),:),1)>tol;
    k_relatives=sum(gram_k,1)>1;
    for i=1:nproj-1
        for j=i+1:nproj            
            ij_map(i,j)=ij_map(i,j)+k_relatives(i)*k_relatives(j);
        end
    end
end
ij_map=triu(ij_map<1);
ij_offset=zeros(nchoosek(nproj,2),1);
prev_offset=0;
idx=0;
for i=1:nproj-1
    for j=i+1:nproj
        idx=idx+1;
        prev_offset=prev_offset+ij_map(i,j);
        ij_offset(idx)=ij_map(i,j)*prev_offset;
    end
end
    








