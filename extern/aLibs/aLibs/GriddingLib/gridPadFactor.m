function pf=gridPadfactor(mode);
% Interprets the string mode and returns the padfactor used in the gridding
% algorithms.
% mode is taken to be either 'grid' --> pf=1.25
% or else 'sinc' --> pf=1.5;
if lower(mode(1))=='g'
    pf=1.25;
elseif lower(mode(1))=='s'
    pf=2;
else
    error(['Invalid mode: ' mode]);
end;
