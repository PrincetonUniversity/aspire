function [ P2,R,R2,x ] = cryo_epsdS_outofcore( instackname,samples_idx,max_d, verbose )
%
% CRYO_EPSDS_OUTOFCORE Estimate the 2D isotropic power spectrum.
%
% [ P2,R,R2,x ] = cryo_epsdS( instackname,samples_idx,max_d, verbose )
%   Estimate the 2D isotropic power spectrum of the projections in the MRC
%   file instackname  using in each image only the pixels given in
%   samples_idx. Typically, samples_idx will correspond to the pixles in
%   the image that are outside a certain raidus (where there is no
%   particle). 
%
% See cryo_epsdS for more details.
%
% Yoel Shkolnisky, May 2016.


if ~exist('verbose','var')
    verbose=0;
end

instack=imagestackReader(instackname);
p=instack.dim(1);
if instack.dim(2)~=p
    error('Images must be square');
end

if max_d>=p
    warning('max_d too large. Setting max_d to %d',p-1);
    max_d=p-1;
end

% Estimate the 1D isotropic autocorrelation function
[R,x,~]=cryo_epsdR_outofcore(instackname,samples_idx,max_d,verbose);

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
for j=1:instack.dim(3)
    im=instack.getImage(j);
    E = E+ sum((im(samples_idx) - mean(im(samples_idx))).^2);
end
meanE=E./(numel(samples_idx)*instack.dim(3));  % Mean energy of the noise samples
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
