
function [arr_tpd]=multi_transpose(arr)
    dims=size(arr);
    if dims(1)~=dims(2)
        disp('2 first dims are not equal ==> matrices are not square');
        arr_tpd=arr;
        return;
    end
    
    l=length(dims);
    if l==2
        d3=[];
    else
        d3=3:l;
    end
    arr_tpd=permute(arr,[2,1,d3]);
end