function rotations=genNearRotations(R,axisdeviation,Naxes,gammadeviation,Ngamma)

% Create rotations near the estimated rotation. Find a plane perpendicular
% to the estimated rotation axis, find points close to the rotation axis on
% this plane, and project them back to the spehere. These are the axes of
% the new rotations to try.

[rotaxis,gamma]=rot_to_axisangle(R);

[Q,~]=qr([rotaxis rand(3,2)]);
assert(norm(Q*Q.'-eye(3))<1.0e-14,'Q not orthogonal'); 

if abs(det(Q)-1)>0.1 % Vector are left handed so flip
    Q(:,[3,2])=Q(:,[2,3]);
end

% Basis for the plane perpendicular to the projection direction rotaxis.
u_x=Q(:,2);
u_y=Q(:,3);

%Naxes=50;  % Generates Naxes new rotation axes
alpha=(rand(2,Naxes)-1/2)*2; % Random coordinates in [-1,1] of points on 
                             % the plane perpendicular to rotaxis.
alpha(:,1)=[0 0].'; % Make the current best guess one of the points
% 0.1 corresponds roughly to a deviation of +/- 5 degrees from the rotation
% axis
newaxes=repmat(rotaxis,1,Naxes)...
    +sin(axisdeviation).*repmat(u_x,1,Naxes).*repmat(alpha(1,:),3,1)...
    +sin(axisdeviation).*repmat(u_y,1,Naxes).*repmat(alpha(2,:),3,1);
newaxes=newaxes./repmat(sqrt(sum(newaxes.^2,1)),3,1); % Project the axes back to the sphere.

% Generate new in-planes rotations that deviate by +/- degrees from the
% estimated in-plane rotation
gamma_degress=gamma/pi*180;
if mod(Ngamma,2)==0
    Ngamma=Ngamma+1; % Make sure to have an odd number of points so that 
                     % the current best guess is one of the points.
end
newgammas=(linspace(-gammadeviation,gammadeviation,Ngamma)+gamma_degress)/180*pi;

rotations=zeros(3,3,Naxes*numel(newgammas));
idx=1;
for kk=1:Naxes
    for jj=1:numel(newgammas)
        rotations(:,:,idx)=axisangle_to_rot(newaxes(:,kk),newgammas(jj));
        idx=idx+1;
    end
end

