function npf_out = gaussian_filter_imgs(npf)
%
% Applies a Guassian filter to the images
%
% Input parameters:
%   npf      A 3D array where each image npf(:,:,i) corresponds to the Fourier
%            transform of projection i.
%
% Output parameters: 
%   npf_out The three-dim array where each slice is filtered

[n_r,n_theta,nImages] = size(npf);
% pf is of size n_rxn_theta. Convert pf into an array of size
% (2xn_r-1)xn_theta, that is, take then entire ray through the origin.
% note that the minus one in 2xn_r-1 is because the dc term appers both in
% an angle theta and in angle theta+pi.
pf=[flipdim(npf(2:end,n_theta/2+1:end,:),1) ; npf(:,1:n_theta/2,:) ];
temp = pf;
% allocate pf to hold a double-cover of all n_theta rays
pf = zeros(2*n_r-1, n_theta, nImages);
pf(:,1: n_theta/2,:)     = flipdim(temp, 1); % the second-cover of all rays are simply heading the opposite side
pf(:,n_theta/2+1:end,:)  = temp;

pf = reshape(pf,2*n_r-1,n_theta*nImages);

rmax = n_r-1; %the dc term appears only once (i.e., in one of the half-rays)
rk = -rmax:rmax; rk=rk(:);
H =  sqrt(abs(rk)).*exp(-rk.^2/(2*(rmax/4).^2));
pf = bsxfun(@times,pf,H);
% the volume is real, hence fourier-trnsform is conjugate-symmetric. Hence
% all we need are half-rays.
pf = pf(1:n_r,:);
pf = flipdim(pf,1); % flip back again by ascending frequency
npf_out = reshape(pf,n_r,n_theta,nImages);

end
