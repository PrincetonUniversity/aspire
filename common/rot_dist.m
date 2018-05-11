function [ dist ] = rot_dist( rot1, rot2)
%
% ROT_DIST  Compute differnce between rotations.
%
% [ dist ] = rot_dist( rot1, rot2)
%   Compute the distance between two sets of rotations.
%   The used norm is the angle of rotation of rot1*rot2^(-1). rot1 and rot2
%   can be either two arrays of quaternions (4 x n1) or two arrays of
%   rotation matrices (3x3xn1). Returns a vector of legnth n1.
%
% Indentical to the function dist_between_rot, but compares two vectors of
% corresponding rotations, as opposed to all rotations from one list with
% all rotations from the other list.
%
% Yoel Shkolnisky, April 2017.
% Based on code by Y. Aizenbud,  5 Jul 2015

%case of quaternions
if size(rot1,1) == 4
    r1 = zeros(3,3,size(rot1,2));
    for ind = 1: size(rot1,2)
        r1(:,:,ind) = q_to_rot(rot1(:,ind));
    end
    r2 = zeros(3,3,size(rot2,2));
    for ind = 1: size(rot2,2)
        r2(:,:,ind) = q_to_rot(rot2(:,ind));
    end
    rot1 = r1;
    rot2 = r2;
end
dist = zeros(size(rot1,3),1);
for ind = 1:size(rot1,3)    
    R = rot1(:,:,ind)*rot2(:,:,ind)^-1;
    %calculating the angle of the rotation. equivalent to a =
    %vrrotmat2vec; dist(ind1,ind2) = a(4)
    mtrc = sum(diag(R));
    dist(ind)= acos((mtrc - 1)/2);      
end

