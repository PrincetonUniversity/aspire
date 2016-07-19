function [ P2,R,R2,x ] = cryo_epsdS( imstack,samples_idx,max_d, verbose )
%
% CRYO_EPSDS Estimate the 2D isotropic power spectrum.
%
% [ P2,R,R2,x ] = cryo_epsdS( imstack,samples_idx,max_d, verbose )
%   Estimate the 2D isotropic power spectrum of a given image stack imstack
%   using in each image only the pixels given in samples_idx. Typically,
%   samples_idx will correspond to the pixles in the image that are outside
%   a certain raidus (where there is no particle).
%   The power spectrum is esimated using the correlogram method.
%
% Input parameters:
%   imstack      Stack of projections. imstack(:,:,k) is the k'th
%                projection in the stack. Projections must be square and
%                can have odd or even dimensions.
%   samples_idx  List of pixel indices to use for autocorrelation
%                estimation.
%   max_d        The autocorrelation is estimated for distances up to
%                max_d. The value of max_d should be much less than
%                the number of data samples, N. A typical value in the
%                literature is N/1, but this may result in a spectrum that
%                is too smooth.
%   verbose      Nonzero to print verbose progress messages. Default 0.
%
% Output parameters:
%   P2           2D power spectrum function. If each image is of size
%                pxp, then P2 is of size (2p-1)x(2p-1). P2 is always real.
%   R            1D isotropic autocorrelation function.
%   R2           2D isotropic autocorrelation function.
%   x            Distances at which the autocorrelction R was estimated.
%
% Yoel Shkolnisky and Mor Cohen, May 2011.
%
% Revised October 2014 (Y.S.):
% 1. Remvoe norm_flag. The power spectrum is always normalized such that
%    its energy is equal to the average energy of the noise samples used
%    to estimate it.
% 2. Remove window_type. Always use bartlett window, as it has positive
%    Fourier transform.
% 3. Remove bias_flag. Not needed.
%
% Revised Y.S. March 3, 2016    Add R2 as output variable.


if ~exist('verbose','var')
    verbose=0;
end

p=size(imstack,1);
if size(imstack,2)~=p
    error('Images must be square');
end

if max_d>=p
    warning('max_d too large. Setting max_d to %d',p-1);
    max_d=p-1;
end

% Estimate the 1D isotropic autocorrelation function
[R,x,~]=cryo_epsdR(imstack,samples_idx,max_d,verbose);

% Use the 1D autocorrelation estimted above to populate an array of the 2D
% isotropic autocorrelction. This autocorrelation is later Fourier
% transformed to get the power spectrum.
R2=zeros(2*p-1); % Two dimensional isotropic autocorrelation
dsquare=x.*x;
for i=-max_d:max_d
    for j=-max_d:max_d
        d=i*i+j*j;
        if d<=max_d*max_d
            idx=bsearch(dsquare,d*(1-1.0e-13),d*(1+1.0e-13));
            % Just subtracting/adding 1.0e-13 from d won't work since for
            % d>1000 this addition/subtraction disappears due to loss of
            % significant digits.
            if isempty(idx) || numel(idx)>1
                error('Something went wrong');
            end
            R2(i+p,j+p)=R(idx);
        end
    end
end


% Window te 2D autocorrelation and Fourier transform it to get the power
% spectrum. Always use the Gaussian window, as it has positive Fourier
% transform.  
w=gwindow(p,max_d);
P2 = cfft2(R2.*w);
err=norm(imag(P2(:)))/norm(P2(:));
if err>1.0e-12
    warning('Large imaginary components in P2 = %d',err);
end
P2=real(P2);

% Normalize the power spectrum P2. The power spectrum is normalized such
% that its energy is equal to the average energy of the noise samples used
% to estimate it.
E=0;  % Total energy of the noise samples used to estimate the power spectrum.
for j=1:size(imstack,3)
    im=imstack(:,:,j);
    E = E+ sum((im(samples_idx) - mean(im(samples_idx))).^2);
end
meanE=E./(numel(samples_idx)*size(imstack,3));  % Mean energy of the noise samples
P2=(P2./sum(P2(:))).*meanE.*numel(P2); 
    % Normalize P2 such that its mean energy is preserved and is equal to
    % meanE, that is, mean(P2(:))==meanE. That way the mean energy does not
    % go down if the number of pixels is artifically changed (say be
    % upsampling, downsampling, or cropping). Note that P2 is already in
    % units of energy, and so the total energy is given by sum(P2(:)) and
    % not by norm(P2(:)).

% Check that P2 has no negative values.
% Due to the truncation of the Gaussian window, we get small negative
% values. So unless they are very big, we just ignore them.
negidx=find(P2<0);
if numel(negidx)>0
    maxnegerr=max(abs(P2(negidx)));
    if verbose
        fprintf('maximal negative value = %e\n',maxnegerr);
    end
    if maxnegerr>1.0e-2
        negnorm=norm(P2(negidx));
        warning('Power specrtum P2 has negative values with energy %e',negnorm);
    end
    P2(negidx)=0; % Zero small negative estimates.
end
