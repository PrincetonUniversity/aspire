function [phi2,mult90]=adjustrotate(phi)
%
% Decompose a rotation CCW by phi into a rotation of mult90 times 90
% degrees followed by a rotation by phi2, where phi2 is between -45 and 45.
% mult90 is an integer between 0 and 3 describing by how many multiples of
% 90 degrees the image should be rotated so that an additional rotation by
% phi2 is equivalent to rotation by phi.
%
% Yoel Shkolnisky, February 2011.

phi=mod(phi,360);

mult90=0;
phi2=phi;

% Note that any two consecutive cases can be combine, but I decided to
% leave them separated for clarity.
if     phi>=45  && phi<90
    mult90=1;   
    phi2=-(90-phi);
elseif phi>=90  && phi<135    
    mult90=1;       
    phi2=phi-90;
elseif phi>=135 && phi<180
    mult90=2;    
    phi2=-(180-phi);
elseif phi>=180 && phi<225
    mult90=2;    
    phi2=phi-180;
elseif phi>=215 && phi<270
    mult90=3;
    phi2=-(270-phi);
elseif phi>=270 && phi<315
    mult90=3;    
    phi2=phi-270;
elseif phi>=315 && phi<360
    mult90=0;
    phi2=phi-360;
end