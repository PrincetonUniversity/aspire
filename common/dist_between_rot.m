function [ dist ] = dist_between_rot( rot1, rot2)
%dist_between_rot(rot1, rot2) returns the distance between two rotations
%   The norm that is used is the angle of rotation of rot1*rot2^(-1)
%   rot1 and rot2 can be either array of quaternions (4 x n1 and 4 x n2) or
%   rotation matrices ( 3x3xn1 and 3x3xn2). the returned value dist is a
%   matrix of size n1 x n2
%
% Written by Y. Aizenbud,  5 Jul 2015
    
    %case of quaternions
    if size(rot1,1) == 4 
        r1 = zeros(3,3,size(rot1,2));
        for ind1 = 1: size(rot1,2)
            r1(:,:,ind1) = q_to_rot(rot1(:,ind1));
        end
        r2 = zeros(3,3,size(rot2,2));
        for ind1 = 1: size(rot2,2)
            r2(:,:,ind1) = q_to_rot(rot2(:,ind1));
        end
        rot1 = r1;
        rot2 = r2;
    end
    dist = zeros(size(rot1,3),size(rot2,3));
    for ind1 = 1:size(rot1,3)
        for ind2 = 1:size(rot2,3)
            R = rot1(:,:,ind1)*rot2(:,:,ind2)^-1;
            %calculating the angle of the rotation. equivalent to a =
            %vrrotmat2vec; dist(ind1,ind2) = a(4)
            mtrc = sum(diag(R));
            dist(ind1,ind2)= acos((mtrc - 1)/2);
                         
        end
    end
end
