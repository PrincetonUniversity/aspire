function ms=QuickShift(m,shift)
% Faster replacement for circshift for time-critical applications.
% About twice as fast as circshift.
% Calls the mex routine QuickShiftMex.
%
sz=size(shift);
if numel(sz)<3 % shift is a 1D or 2D vector: use quick shift code
    if numel(shift)<2
        shift(2)=0;
    end;
    dims=size(m);
    shift=mod(shift,dims);
    ms=QuickShiftMex(m,shift);
else
    ms=circshift(m,shift);
end;
    