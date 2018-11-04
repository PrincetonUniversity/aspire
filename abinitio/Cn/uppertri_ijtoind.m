function ind = uppertri_ijtoind(i,j,n)
% (1,2) ->1, (n-1,n)->n*(n-1)/2

assert(i<j,'i must be smaller than j');

ind = (2*n-i)*(i-1)/2+j-i;

end