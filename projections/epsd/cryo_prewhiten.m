function [ proj, filter, nzidx] = cryo_prewhiten(proj, noise_response, rel_threshold)
% Pre-whiten a stack of projections using the power spectrum of the noise.
%  
%  Input:
%   proj            Stack of images.
%   noise_response  2d image with the power spectrum of the noise. If all
%      images are to be whitened with respect to the same power spectrum,
%      this is a single image. If each image is to be whitened with respect
%      to a different power spectrum, this is a three-dimensional array with
%      the same number of 2d slices as the stack of images.
%   rel_threshold   The relative threshold used to determine which frequencies
%                   to whiten and which to set to zero. If empty (the default)
%                   all filter values less than 100*eps(class(proj)) are
%                   zeroed out, while otherwise, all filter values less than
%                   threshold times the maximum filter value for each filter
%                   is set to zero.
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

if nargin < 3
    rel_threshold = [];
end

delta=eps(class(proj));

n=size(proj, 3);
L=size(proj, 1);
l=floor(L/2);
K=size(noise_response, 1);
k=ceil(K/2);

if size(noise_response, 3) ~= n && size(noise_response, 3) ~= 1
    error('The number of filters must be either 1 or n.');
end

% The whitening filter is the sqrt of of the power spectrum of the noise.
% Also, normalized the enetgy of the filter to one.
filter=sqrt(noise_response);      
filter=filter/norm(filter(:));

% The power spectrum of the noise must be positive, and then, the values
% in filter are all real. If they are not, this means that noise_response
% had negative values so abort.
assert(norm(imag(filter(:)))<10*delta*size(filter, 3)); % Allow loosing one digit.
filter=real(filter);  % Get rid of tiny imaginary components, if any.

% The filter should be cicularly symmetric. In particular, it is up-down
% and left-right symmetric.
assert(norm(filter-flipud(filter))<10*delta); 
assert(norm(filter-fliplr(filter))<10*delta);

% Get rid of any tiny asymmetries in the filter.
filter=(filter+flipud(filter))./2;
filter=(filter+fliplr(filter))./2;

% The filter may have very small values or even zeros. We don't want to
% process these so make a list of all large entries.
if isempty(rel_threshold)
    nzidx = find(filter>100*delta);
else
    nzidx = find(bsxfun(@gt, filter, rel_threshold*max(max(filter, [], 1), [], 2)));
end

nzidx = nzidx(:);

fnz=filter(nzidx);

% Pad the input projections
pp=zeros(K);
p2=zeros(L,L,n);
for idx=1:n
    if mod(L,2)==1 % Odd-sized image
        pp(k-l:k+l, k-l:k+l)=proj(:,:,idx); % Zero pad the image to twice the size.
    else
        pp(k-l:k+l-1, k-l:k+l-1)=proj(:,:,idx); % Zero pad the image to twice the size.
    end
    
    fp=cfft2(pp); % Take the Fourier transform of the padded image.
    p=zeros(size(fp));

    % Divide the image by the whitening filter, 
    % but onlyin places where the filter is
    % large. In frequnecies where the filter is
    % tiny  we cannot pre-whiten so we just put
    % zero.
    p(nzidx) = bsxfun(@times, fp(nzidx), 1./fnz);
    pp2 = icfft2(p); % pp2 for padded p2.
    assert(norm(imag(pp2(:)))/norm(pp2(:))<1.0e-13); % The resulting image should be real.
    
    if mod(L,2)==1
        p2(:,:,idx) = pp2(k-l:k+l, k-l:k+l);
    else
        p2(:,:,idx) = pp2(k-l:k+l-1, k-l:k+l-1);
    end

end
proj = real(p2);

