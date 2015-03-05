function printProgressBarHeader

global ticsprinted
ticsprinted=0;

Ndots=85;
textbar=['0%' repmat('.',1,floor(Ndots/4)) '25%' ...
    repmat('.',1,floor(Ndots/4)) '50%'...
    repmat('.',1,floor(Ndots/4)) '75%'... 
    repmat('.',1,Ndots-3*floor(Ndots/4)) '100%'];
fprintf('%s\n\n',textbar);