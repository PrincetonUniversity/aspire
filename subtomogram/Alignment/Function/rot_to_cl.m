function [ ind ] = rot_to_cl(rot_est, tiltq, n_theta, ns)
% This function takes a rotation between two subtomograms and compute the
% common line indices on the subtomogram slices

% Input:
% rot_est: the rotation between two subtomograms in rotation matrix 
% tiltq: 4*(# of tilts in one subtomogram) is a matrix representing the
% tilt series tilted around x-axis with z axis as the projection angle,
% where each column is the quaternion for the rotation of one of the
% slices.
% n_theta: the number of angles that 2pi is divided into, i.e., the
% resolution of the angles. Try doubling the original grid number of
% rotations.
%
% Output: 
% ind: the index matrices indicating which radial line is the common line
% in the two subtomograms. ind(i,j) is the index of the common line on the
% i-th slice from the j-th slice, assumming that the index 1 radial line is
% the common line shared by all slices within one subtomogram. The first
% subtomograms is the one we use as the reference coordination, and the
% slice number ranges from 1 to 9, starting from the closest to the 
% positive x-y plane. While the slices in second subtomogram has index 9-18
% in the matrix.

%ns = 9; % number of tilt series in a subtomogram
ind = ones(2*ns) - eye(2*ns);
v0 = [0 0 1]'; % z axis unit vector represent the vector orthogonal to the 
                % original xy plane
dis_theta = 2*pi/n_theta;
rot_est = rot_est';

% compute the column vectors orthogonal to the first and second subtomogram 
% slices
V1 = zeros(3,ns); 
V2 = zeros(3,ns);
for i = 1:ns
    rot = q_to_rot(tiltq(:,i));
    rot = rot';
    V1(:,i) = rot*v0;
    V2(:,i) = rot_est*rot*v0;
end

% compute their intersection lines
% discretize the xy-plane unit circle

% rotate the original unit circle to generate the tilt series of the two
% subtomograms and find their intersecting lines
% finding their intersecting line by calculating the point along the
% circular cone and see where the two surfaces have the same z value
nx = [-1 0 0]'; % the orientation is looking from the x axis 
for i = 1:ns % index of the first subtomogram
    for j = 1:ns % index of the second subtomogram
        n = cross(V1(:,i),V2(:,j));
        n = n/norm(n); % common line
        n2 = rot_est*nx; % starting line on the slice
        mid = (ns+1)/2; 
        % viewing direction from the -z-axis, the orientation of the line
        % numbers on each slice changed from clockwise to counter-clockwise
        theta1 = acos(dot(n,nx))*sign(dot(cross(n,nx),-V1(:,mid))); 
        theta2 = acos(dot(n,n2))*sign(dot(cross(n,n2),-V2(:,mid)));
        if theta1 < 0
            ind(i,j+ns) = n_theta + 1 + round(theta1/dis_theta);
            if round(theta1/dis_theta) == 0 %in order to avoid 73
                ind(i,j+ns) = 1;
            end
        else
            ind(i,j+ns) = 1 + round(theta1/dis_theta);
        end
        if theta2 < 0
            ind(j+ns,i) = n_theta + 1 + round(theta2/dis_theta);
            if round(theta2/dis_theta) == 0 %in order to avoid 73
                ind(j+ns,i) = 1;
            end
        else
            ind(j+ns,i) = 1 + round(theta2/dis_theta);
        end
    end
end
end

% Test with rot_s = the rotation of z->-y

%     0     0    -1
%     0     1     0
%     1     0     0

