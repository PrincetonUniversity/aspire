function [noisy_cl,mask]=perturb_clmatrix(clmatrix,L,p)
%
% Perturb the common lines matrix
%
% clmatrix is a clean common lines matrix that was constructed using
% angular resolution of L radial lines per image.
% p is the probability that correlating two projection images yields the
% correct common line.
%
% noisy_cl is the perturb common lines matrix. mask is a 0-1 matrix of the
% same size as the common lines matrix, where mask(k1,k2) with k1<k2 is 1
% if the common line between images k1 and k2 was perturbed randomly, and 0
% otherwise. The lower triangular part of mask is zero.

K=size(clmatrix,1);
mask = zeros(K); % to verify the success of the algorithm we store which common lines were perturbed
noisy_cl=clmatrix;
for k1=1:K;
    for k2=(k1+1):K;
        r = rand;
        % with probability 1-p the common line needs to be perturbed
        if (r > p)
            mask(k1,k2) = 1;
            noisy_cl(k1,k2) = floor(rand*L)+1;
            noisy_cl(k2,k1) = floor(rand*L)+1;
        end;
    end;
end;