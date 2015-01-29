function [rotaxis,gamma]=rot_to_axisangle(R)
%
% ROT_TO_AXISANGLE  Convert rotation matrix to axis-angle representation
%
% [rotaxis,gamma]=rot_to_axisangle(R)
%    Represent the rotation matrix R as a rotation by angle gamma around
%    the axis rotaxis
%
% Yoel Shkolnisky, January 2015

cos_gamma=(trace(R)-1)/2;

if abs(cos_gamma-1)<1.0e-8  % in-plane rotation is near zero so there is
    % no rotation at all.
    gamma=0;
    rotaxis=[0;0;1]; % Can return any vector, but make it unit norm for the 
                     % assertion below.
elseif abs(cos_gamma+1)<1.0e-8 % cos_gamma_ref is near -1
    gamma=pi;
    B=(eye(3)+R)/2;
    Bvec=[B(1,1);B(2,2);B(3,3)];
    if sum(Bvec>0)==3 % All Bii>0
        rotaxis=sqrt(Bvec);
    else
        Bvec=[B(1,2);B(1,3);B(2,3)];
        if sum(Bvec>0)==1 && sum(Bvec<0)==2
            rotaxis=[-sqrt(B(1,1));+sqrt(B(2,2));+sqrt(B(3,3))];
        elseif sum(abs(Bvec)<1.0e-8)==2
            rotaxis=[0;+sqrt(B(2,2));+sqrt(B(3,3))];
            % The following is also correct in this case
            % rotaxis_ref=[-sqrt(B(1,1));+sqrt(B(2,2));-sqrt(B(3,3))];
        elseif sum(abs(Bvec)<1.0e-8)==3
            rotaxis=[sqrt(B(1,1));0;0];
        end
    end
else
    gamma=acos(cos_gamma);
    N=(R-R.')./(2*sin(gamma));
    assert(norm(N+N.')<1.0e-13) % Should be skew symmetric
    rotaxis=[N(3,2);N(1,3);N(2,1)];
end

assert(abs(norm(rotaxis)-1)<1.0e-8);

% Just to make sure, reconstruct the rotation matrix from the axis and the
% angle and compare to the input rotation.
R2=axisangle_to_rot(rotaxis,gamma);    
assert(norm(R-R2)<1.0e-13,...
    'Failed to extract parameters of rotation');