function [ O1, O2 ] = bispec_Operator(freq)
%   It generates the operators for bispectrum computation
%   Input: 
%       freq: frequencies computed from FBsPCA.
%   Output: 
%       O1 and O2 matrices used for computing bispectrum
%Zhizhen Zhao 09/01/2012

max_freq=max(freq);
ind=1;
list=zeros(100000, 3);
for i1=1:max_freq
        for j1=1:i1-1
                k1=i1+j1;
                id1=find(freq==i1);
                id2=find(freq==j1);
                id3=find(freq==k1);
                if ~isempty(id3)
                    for k1=1:length(id1)
                        for k2=1:length(id2)
                            for k3=1:length(id3)
                                list(ind, :)=[id1(k1), id2(k2), id3(k3)];
                                ind=ind+1;
                            end;
                        end;
                    end;
                end;
        end;
end;

list=list(1:ind-1, :);
val=ones(size(list, 1), 3);
val(:, 3)=-1;
n_col=size(list, 1); 
list=list(:);
col=repmat([1:n_col]', 3, 1);
O1=sparse(col, list,ones(length(list), 1), n_col, length(freq));
O2=sparse(col, list,val(:), n_col, length(freq));

end

