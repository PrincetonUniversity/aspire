function [rotaxis,gamma]=rot_to_axisangle(R)
%
% ROT_TO_AXISANGLE  Convert rotation matrix to axis-angle representation
%
% [rotaxis,gamma]=rot_to_axisangle(R)
%    Represent the rotation matrix R as a rotation by angle gamma around
%    the axis rotaxis
%
% Yoel Shkolnisky, January 2015
%
% Revisions
% Y.S. May 17, 2015     Fix the case of gamma=pi. Test which component is
%      zero using a tolerance delta instead of 0. Otherwise, if some
%      component is 10^-16 it will be conisdered non-zero resulting in a
%      too large error in R2 compared to R. 

cos_gamma=(trace(R)-1)/2;

if abs(cos_gamma-1)<1.0e-8  % in-plane rotation is near zero so there is
    % no rotation at all.
    gamma=0;
    rotaxis=[0;0;1]; % Can return any vector, but make it unit norm for the 
                     % assertion below.
elseif abs(cos_gamma+1)<1.0e-8 % cos_gamma_ref is near -1
    gamma=pi;
    B=(eye(3)+R)/2;
    Bvec=[B(1,2);B(2,3);B(1,3)];
    delta=1.0e-8; % A entry b in B is considered zero is abs(b)<delta.
    if sum(Bvec>delta)==3 % All Bii>0
        rotaxis=sqrt(diag(B));
    elseif sum(Bvec>delta)==1 && sum(Bvec<-delta)==2
            if      B(1,2)<-delta && B(1,3)<-delta % 1 is the common index
                rotaxis=[-sqrt(B(1,1));+sqrt(B(2,2));+sqrt(B(3,3))];
            elseif  B(1,2)<-delta && B(2,3)<-delta % 2 is the common index
                rotaxis=[+sqrt(B(1,1));-sqrt(B(2,2));+sqrt(B(3,3))];
            elseif  B(1,3)<-delta && B(2,3)<-delta % 3 is the common index
                rotaxis=[+sqrt(B(1,1));+sqrt(B(2,2));-sqrt(B(3,3))];
            else
                error('No common index found - should never happen');
            end            
    elseif sum(abs(Bvec)<delta)==2
        % In this case exactly one of B(1,1), B(2,2), B(3,3) is zero.
            if      abs(B(1,1))<delta
                rotaxis=[0;+sqrt(B(2,2));+sqrt(B(3,3))];
                rotaxis(3)=rotaxis(3)*sign(B(2,3));
            elseif  abs(B(2,2))<delta
                rotaxis=[+sqrt(B(1,1));0;+sqrt(B(3,3))];
                rotaxis(3)=rotaxis(3)*sign(B(1,3));
            elseif  abs(B(3,3))<delta
                rotaxis=[+sqrt(B(1,1));+sqrt(B(2,2));0];
                rotaxis(2)=rotaxis(2)*sign(B(1,2));
            else
                error('Exactly one of B(1,1),B(2,2), B(3,3) should have been zero?!');
            end
    elseif sum(abs(Bvec)<delta)==3
        if      abs(B(1,1))>delta
            rotaxis=[sqrt(B(1,1));0;0];
        elseif  abs(B(2,2))>delta
            rotaxis=[sqrt(B(2,2));0;0];
        elseif  abs(B(3,3))>1.0e-8
            rotaxis=[sqrt(B(3,3));0;0];
        else
            error('Exactly one of B(1,1),B(2,2), B(3,3) should have been NON-zero?!');
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
if norm(R-R2)>=1.0e-8
    aaa=1;
end

assert(norm(R-R2)<1.0e-8,...
    'Failed to extract parameters of rotation');