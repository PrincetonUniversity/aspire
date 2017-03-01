function clerr=syncconsistency(rotations,clmatrix,L)
%
% Wrapper for backward compitability. 
%
% Yoel Shkolnisky, March 2015

warning('This function is depreacted. Call cryo_syncconsistency instead.');
clerr=cryo_syncconsistency(rotations,clmatrix,L);
