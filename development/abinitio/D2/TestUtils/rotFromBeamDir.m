
function [R]=rotFromBeamDir(v)

R=zeros(3,3);
v=v/norm(v);
R(:,3)=v;
R(:,2)=[-v(2);v(1);0]/norm([-v(2);v(1);0]);
R(:,1)=cross(R(:,2),R(:,3))/norm(cross(R(:,2),R(:,3)));
end