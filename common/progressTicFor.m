function progressTicFor(j,N)
%
% DOES NOT WORK FOR PARFOR.

global ticsprinted

if isempty(ticsprinted)
    ticsprinted=0;
end

fractioncompleted=floor(j*100/N);

if fractioncompleted>ticsprinted
    for k=1:fractioncompleted-ticsprinted
        fprintf('\b|\n');
    end
    ticsprinted=fractioncompleted;
end
