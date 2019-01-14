
function [idx_map]=linIdx2Trip(n)
    idx_map=zeros(nchoosek(n,3),3);
    idx=0;
    for i=1:n-2
        for j=i+1:n-1
            for k=j+1:n
                idx=idx+1;
                idx_map(idx,:)=[i,j,k];
            end
        end
    end
end