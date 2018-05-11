function [ O1, O2 ] = bispec_Operator_1(freq)
%Description:
%   It generates the operators for bispectrum computation
%   Input: 
%       freq: frequencies computed from FBsPCA.
%   Output: 
%       O1 and O2 matrices used for computing bispectrum. O1 for the amplitude and O2 for the phase. 
%Zhizhen Zhao 09/01/2012
%  Updated Jan 2015 for complex imag

max_freq=max(freq);
count = 0;
list=zeros(100000, 3);
for i1=2:max_freq-1
        for j1=1:min(i1-1, max_freq-i1)
                k1=i1+j1;
                id1=find(freq==i1);
                id2=find(freq==j1);
                id3=find(freq==k1);
                nd1 = length(id1);
                nd2 = length(id2);
                nd3 = length(id3);
                nd = nd1*nd2*nd3;
                if nd~=0
                    tmp1 = repmat(id1, nd2, 1);
                    tmp2 = kron(id2, ones(nd1, 1));
                    tmp = [ tmp1, tmp2 ];
                    tmp1 = repmat(tmp, nd3, 1);
                    tmp3 = kron(id3, ones(nd1*nd2, 1));
                    list( count+1 : count + nd, : ) = [tmp1, tmp3];
                    count = count + nd;
                end;
        end;
end;

list=list(1:count, :);
val = ones(size(list));
val(:, 3) = -1;
n_col=size(list, 1); 
list=list(:);
col=repmat([1:n_col]', 3, 1);
O1 = sparse(col, list, ones(length(list), 1), n_col, length(freq));
O2 = sparse(col, list, val, n_col, length(freq));   
end

