
function [idx_map]=lin2sub3_map(N)
idx_map=zeros(nchoosek(N,3),3);
idx=0;
for i=1:N-2
    for j=i+1:N-1
        for k=j+1:N
            idx=idx+1;
            idx_map(idx,:)=[i,j,k];
        end
    end
end