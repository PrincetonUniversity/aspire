function pf2=cryo_raynormalize(pf)
%
% Normalize a dataset of Fourier rays so that each ray has energy 1.
% 
% pf is a 3D array of the Fourier transform of the projections. 
% pf(:,:,k) is the polar Fouier transform of the k'th projection.
%

n_proj=1;
if ndims(pf)==3
    n_proj=size(pf,3);
end

n_theta=size(pf,2);

pf2=pf; % create a copy of the data for normalization purposes
for k=1:n_proj
    for j=1:n_theta
        nr=norm(pf2(:,j,k));
        if nr<1.0e-13
            warning('Ray norm is close to zero. k=%d  j=%d',k,j);
        end
        pf2(:,j,k)=pf2(:,j,k)./nr;
    end
end