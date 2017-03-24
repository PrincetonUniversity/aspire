function org=fctr(n)
% function org=fctr(n)
% Center of an FFT-shifted image.  We use this center coordinate for all
% rotations and centering operations.
% If n is even, org=[n/2+1;n/2+1].
% If n is odd, org is [(n+1)/2;(n+1)/2].
% n can be a two-element vector n=[nx ny]
if numel(n)<2
    n=[n n];
end;
org=ceil((n+1)/2);

