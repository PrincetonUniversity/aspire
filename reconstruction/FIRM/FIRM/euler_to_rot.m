function [ R ] = euler_to_rot( angle )
% in spider, The rotation matrices used are defined as:
% (http://www.wadsworth.org/spider_doc/spider/docs/euler.html)
% v = Rv', where  R is the matrix for transforming vector v' to vector v.
%
% R = R_psi * R_theta * R_phi
%
%
% R_psi   =   cos_psi  sin_psi  0
%             -sin_psi  cos_psi  0
%                 0         0      1
%
%
% R_theta =   cos_theta  0 -sin_theta
%                 0        1     0
%              sin_theta  0  cos_theta
%
%
% R_phi   =   cos_phi  sin_phi  0
%             -sin_phi  cos_phi  0
%                 0         0      1

angle=angle(:)';
cos_ang=cos(angle/360*2*pi);
sin_ang=sin(angle/360*2*pi);
cos_psi=cos_ang(1);
cos_theta=cos_ang(2);
cos_phi=cos_ang(3);
sin_psi=sin_ang(1);
sin_theta=sin_ang(2);
sin_phi=sin_ang(3);


R_psi   =  [ cos_psi  sin_psi  0;
    -sin_psi  cos_psi  0;
    0         0      1];


R_theta =   [cos_theta  0 -sin_theta;
    0        1     0;
    sin_theta  0  cos_theta];


R_phi   =   [cos_phi  sin_phi  0 ;
    -sin_phi  cos_phi  0 ;
    0         0      1];

R = R_psi * R_theta * R_phi ;
end