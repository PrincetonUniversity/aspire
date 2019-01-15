

% This function computes which pairs of rotations are both equators in D2. 
% Input: 
% eq_filter_angle = which distance from equator is considered an equator
% filterAllEqPairs = '1' if we wish to mark all equator pairs, '0' if we
% mark only equators relative to the same symmetry axis
% addFilter = additional filter table (lower triangular matrix)
%
% Output:
% pairs = Filtered index pairs 
function [pairs,eq_idx,eq_class]=...
    getNonEqPairsIdx(rots,eq_filter_angle,filterAllEqPairs,addFilter)

nrot=size(rots,3);
projDirs=squeeze(rots(:,3,:));
[eq_idx,eq_class,~]=markEquators(projDirs,eq_filter_angle);
eq_idx=single(eq_idx);
eq_table_idx=eq_idx*eq_idx';%places of eq pairs in the table

if filterAllEqPairs
    eq2eq_table=tril(ones(nrot,nrot)-...
        eq_table_idx,-1);
else
    eq2eq_table=tril(ones(nrot,nrot)-...
        eq_table_idx.*(abs(round(eq_class-eq_class'))==0),-1);
end

if nargin>3
   eq2eq_table=eq2eq_table.*addFilter;
end

[pairsI,pairsJ]=meshgrid(1:nrot,1:nrot);
idxI=pairsI(logical((pairsJ>pairsI).*eq2eq_table));
idxJ=pairsJ(logical((pairsJ>pairsI).*eq2eq_table));
pairs=[idxI,idxJ];


