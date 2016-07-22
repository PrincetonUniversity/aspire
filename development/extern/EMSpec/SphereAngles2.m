function angles=SphereAngles2(ntheta,npsi)
% function [angles inds]=SphereAngles(ntheta,npsi)
% Create an array 3 x t of angle values which are ~uniform on a hemisphere,
%   At each theta value,the number
% of computed phi values is approximately nphi*sin(theta)
% 
% % Set up the list of angles
% angles=SphereAngles(18,72);  % 5 degree steps, 813 values
% angles=SphereAngles( 9,36);  % 10 degree steps, 199 values
% 
% Modified so that the theta values go from 0:dtheta:1-dtheta/2 x pi/2.
% fs 10 Mar 2012


nthetaM1=ntheta-1;  % there will be 1 more than this number of theta steps

angles=[];
t=0;
dtheta=pi/2/(nthetaM1+.5);

for i=1:ntheta
    % Assigning theta and psi values:
    theta=(i-1)*dtheta;

    NPsiSteps=min(npsi,sin(theta)*npsi+1);
    for j=1:NPsiSteps
        % psi values are (j-1)*NPsiSteps*2*pi for j=1..NPsiSteps
        psi=(j-1)/NPsiSteps*2*pi;
        % nangles=nangles+1;
        t=t+1;
        angles(t,:)=[0 theta psi];
    end
end;
