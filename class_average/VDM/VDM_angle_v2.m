function [ angle ] = VDM_angle_v2(V, list)
%Compute rotational alignment between data points using eigenvectors from
%VDM
%   Detailed explanation goes here
%   Input:
%       V: 
%           Matrix of size Pxn_eig. P data points, n_eig eigenvectors. Each
%           column is weighted by the eigenvalues.
%       list: 
%           nearest neighbor list.
%   Output: 
%       angle:
%           in-plane rotational alignment bewteen nearest neighbors.
%
%Zhizhen Zhao Aug 2013

angle=sum(V(list(:, 2), :).*conj(V(list(:, 1), :)), 2);
angle=atan2(imag(angle), real(angle));
angle=angle*180/pi;

end

