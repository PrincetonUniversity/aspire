function [ proj ] = Prewhiten_image2d(proj, noise_response )
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

delta=eps(class(proj));

n=size(proj, 3);
L=size(proj, 1); 
l=floor(L/2);
K=size(noise_response, 1);
k=ceil(K/2);

% The whitening filter is the sqrt of of the power spectrum of the noise.
% Also, normalized the enetgy of the filter to one.
filter=sqrt(noise_response);      
filter=filter/norm(filter(:));

% The power spectrum of the noise must be positive, and then, the values
% in filter are all real. If they are not, this means that noise_response
% had negative values so abort.
assert(norm(imag(filter(:)))<10*delta); % Allow loosing one digit.
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
nzidx=find(filter>100*delta);
fnz=filter(nzidx);

fprintf('Whitening...\n');
printProgressBarHeader;
parfor i=1:n
    progressTic(i,n);
    pp=zeros(K);
    if mod(L,2)==1 % Odd-sized image
        pp(k-l:k+l, k-l:k+l)=proj(:, :, i); % Zero pad the image to twice the size.
    else
        pp(k-l:k+l-1, k-l:k+l-1)=proj(:, :, i); % Zero pad the image to twice the size.
    end
    fp=cfft2(pp); % Take the Fourier transform of the padded image.
    p=zeros(size(fp));
    p(nzidx) = fp(nzidx)./fnz; % Divide the image by the whitening filter, 
                               % but onlyin places where the filter is
                               % large. In frequnecies where the filter is
                               % tiny  we cannot pre-whiten so we just put
                               % zero.
    p2 = icfft2(p);
    assert(norm(imag(p2(:)))/norm(p2(:))<1.0e-13); % The resulting image should be real.
    if mod(L,2)==1
        p2 = p2(k-l:k+l, k-l:k+l);   
    else
        p2 = p2(k-l:k+l-1, k-l:k+l-1);
    end
    proj(:, :, i)=real(p2);
end;

end

