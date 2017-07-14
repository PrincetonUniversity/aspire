function [ proj, filter, nzidx] = cryo_prewhiten(proj, noise_response )
% Pre-whiten a stack of projections using the power spectrum of the noise.
%  
%  Input:
%   proj            Stack of images.
%   noise_response  2d image with the power spectrum of the noise.
%
%  Output:
%   proj    Pre-whitened stack of images.
%
% Yoel Shkolnisky and Zhizhen Zhao, July 2013.
% Revisions:
%   08/04/2015  Y.S.   Change all thresholds to support both single and
%                      double precision.
%
%   29/05/2016  Y.S.    Rename Prewhiten_image2d to cryo_prewhiten

delta=eps(class(proj));

n=size(proj, 3);
L=size(proj, 1); 
K=size(noise_response, 1);

% The whitening filter is the sqrt of of the power spectrum of the noise.
% Also, normalized the enetgy of the filter to one.
filter=sqrt(noise_response);      
filter=filter/norm(filter(:));

% The power spectrum of the noise must be positive, and then, the values
% in filter are all real. If they are not, this means that noise_response
% had negative values so abort.
assert(norm(imag(filter(:)))<10*delta); % Allow loosing one digit.
filter=real(filter);  % Get rid of tiny imaginary components, if any.

% The filter should be cicularly symmetric. In particular, should have
% reflection symmetry.
assert(norm(filter-fourier_flip(filter, [1 2]))<10*delta); 

% Get rid of any tiny asymmetries in the filter.
filter = 0.5*(filter + fourier_flip(filter, [1 2]));

% The filter may have very small values or even zeros. We don't want to
% process these so make a list of all large entries.
nzidx=find(filter>100*delta);

nzidx = nzidx(:);

fnz=filter(nzidx);

% TODO: Batch this to avoid memory trouble with large n if the noise_response
% is much larger than the projections.

fprintf('Whitening...\n');
pp = zero_pad(proj, K*ones(1, 2));
fp=cfft2(pp); % Take the Fourier transform of the padded image.
fp = reshape(fp, [L^2 n]);
p=zeros([L^2 n]);
% Divide the image by the whitening filter, 
% but onlyin places where the filter is
% large. In frequnecies where the filter is
% tiny  we cannot pre-whiten so we just put
% zero.
p(nzidx,:) = bsxfun(@times, fp(nzidx,:), 1./fnz);
p = reshape(p, [L*ones(1, 2) n]);
p2 = icfft2(p);
assert(norm(imag(p2(:)))/norm(p2(:))<1.0e-13); % The resulting image should be real.
p2 = extract_center(p2, L*ones(1, 2));
proj = real(p2);

end

