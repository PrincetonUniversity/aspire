
function [idx_map]=lin2sub_map(N)
    idx_map=zeros(nchoosek(N,2),2);
    idx=0;
    for i=1:N-1
        for j=i+1:N
            idx=idx+1;
            idx_map(idx,:)=[i,j];
        end
    end
end