function [estdx,err]=register_translations_2d(im1,im2,refdx,debug,doplot)
%
% REGISTER_TRANSLATIONS_2D  Estimate relative shift between two images.
%
% register_translations_2d(im1,im2,refdx,doplot)
%   Estimate the relative shift between two images im1 and im2 to subpixel
%   accuracy. The function uses phase correlation to estimate the relative
%   shift to within one pixel accuray, and then refines the estimation
%   using Newton iterations.
%   
%   Input parameters:
%   im1,im2 Two images to register. Images must be odd-sized.
%   refidx  Two-dimensional vector with the true shift between the images,
%           used for debugging purposes. (Optional)
%   doplot  If nonzero then the function prints debugging images.
%           (Optional).
%
%   Output parameters
%   estidx  A two-dimensional vector of how much to shift im2 to aligh it
%           with im1. 
%   err     Difference between estidx and refidx.
%
% Yoel Shkolnisy, January 2014.

if ~exist('debug','var');
    debug=0;
end

if exist('refdx','var')
    debug=1;
else 
    refdx=[];    
end

if ~exist('doplot','var')
    doplot=0;
end


% Take Fourer transform of both images and compute the phase correlation
% factors.
hats1 = fftshift( fft2( ifftshift(im1))); % Compute centered Fourier transform
hats2 = fftshift( fft2( ifftshift(im2)));
rhat=hats1.*conj(hats2)./(abs(hats1.*conj(hats2)));
if any(isnan(rhat(:))) || any(isinf(rhat(:)))
    error('phase factors have nans or infs');
end

n=size(im1,1);

% % % if mod(n,2)~=1
% % %     error('Images must be odd-sized');
% % % end
ll=fix(n/2);
freqrng=-ll:n-ll-1;

%doplot=1;
if ~isempty(refdx)
    % Plot true and estimated phase factors    
    [omega_x,omega_y]=ndgrid(freqrng,freqrng);
    omegax=exp(-1i.*2*pi.*omega_x.*refdx(1)/n);
    omegay=exp(-1i.*2*pi.*omega_y.*refdx(2)/n);

    true_phases=omegax.*omegay;
    % The minus is because hats2 is conjugated in rhat , which flips the
    % sign of  its phases.
    
    if doplot
        figure;
        subplot(2,2,1);
        imagesc(real(rhat))
        title('rhat - real part')
        subplot(2,2,2);
        imagesc(real(true_phases));
        title('true phases - real part');
        
        subplot(2,2,3);
        imagesc(imag(rhat))
        title('rhat - imaginary part')
        subplot(2,2,4);
        imagesc(imag(true_phases))
        title('true phases - imaginary part');
    end
    
    % If n is even, then discard the frequnecy which has no conjugate
    % symmetry.
    if mod(n,2)==0
        rhat(:,1)=1;
        rhat(1,:)=1;
        true_phases(:,1)=1;
        true_phases(1,:)=1;
    end
    
    fprintf('L2 error in estimting rhat = %e\n',norm(rhat(:)-true_phases(:))/norm(true_phases(:)));
end


% Compute the relative shift between the images to to within 1 pixel
% accuracy.

% mm is a window function that can be applied to the images before
% computing their relative shift. Experimenting with the code shows that
% windowing does not improve accuracy.

%mm=fuzzymask(n,2,floor(n/4),ceil(n/20));
mm=1; % Windowing does not improve accuracy.
r=fftshift(ifft2(ifftshift(rhat.*mm)));
[~,ii]=max(r(:));  % Find the peak of the IFFT of the phase factors.
cX=fix(n/2)+1; %Find the center
cY=fix(n/2)+1; %Find the center

[sX,sY]=ind2sub(size(r),ii);
est_shift=[cX-sX,cY-sY];

if debug
    fprintf('estimated shift = (%d,%d)\n',est_shift(1),est_shift(2));
end


% Refine shift estimation using Newton iterations
estdx=-1;
MAXITER=200; % Maximal number of Newton iterations
eps=1.0e-8;  % Error to terminate Newton
iter=1;      % Iteration number
x=est_shift; % Initialize Newton from the phase correlation estimated shifts.
x=x(:);
p=1;         % Newton step size

% Use only phases close the origin. These phases are more reliable.
radius=fix(n/2)*0.5;
[xx,yy]=ndgrid(freqrng,freqrng);
idx=find(xx.^2+yy.^2<radius^2); % idx are the indices of the reliable phases.
% Note that we never use the top-most and left-most frequnecies of rhat
% (the frequnecies which have no conjugate-symmetric conuterpart) since in
% such a case, for n even, we would need to modify these values in the
% estimated phases to be real. In other words, for n even, these
% frequnecies do not contain the required information, so ignore them.

if debug && doplot
    % Computed an image containing only the retained phases.
    masked_rhat=zeros(size(rhat));
    masked_rhat(idx)=rhat(idx);
    
    figure;
    subplot(1,2,1);
    imagesc(real(masked_rhat)); axis image;
    subplot(1,2,2);
    imagesc(imag(masked_rhat)); axis image;
end

% A function to evaluate the L2 error of the estimated shifts.
func=@(x,rhat,n,idx) sum(abs(E2(x,rhat(idx),n,idx)).^2);

f=sum(abs(E2(x,rhat(idx),n,idx)).^2);
ferr=abs(f);
while iter<=MAXITER && ferr>eps &&  norm(p)>eps
    df=sum(dE2(x,rhat(idx),n,idx));
    d2f=reshape(sum(d2E2(x,rhat(idx),n,idx),1),2,2);
    % Debug - compute the second derivative numerically by central
    % difference to check the output of d2E2.
    %d2f1=(sum(dE(j,x+1.0e-4,rhat(idx),(n-1)/2,idx))-sum(dE(j,x-1.0e-4,rhat(idx),(n-1)/2,idx)))./(2*1.0e-4)
    
    df=df(:);
    p=-d2f\df;    
    
    % Line search. Instead of taking a full Newton step from x in the
    % direction p, find the minimum of f starting from x in the direction
    % p. Note that using line search is not the most efficient way (see
    % Numerical Receipes).
    fold=f;   
    [xmin,fmin]=linmin(x,p,func,{rhat,n,idx},1.0e-8,0);
    x=xmin;
    f=fmin;
    
%     xold=x;
%     lambda=1;    
%     x=x+p; % Since we have a double root at the minimum.
%     f=sum(abs(E2(x,rhat(idx),n,idx)).^2);
%         
%     while f>fold && lambda>1.0e-4
%         lambda=lambda/2;
%         x=xold+lambda*p;
%         f=sum(abs(E2(x,rhat(idx),n,idx)).^2);
%     end

    ferr=abs(f-fold)/max(abs(f),1);
    iter=iter+1;
end

if debug
    fprintf('On iteration %d, ',iter-1);
end

if iter>=MAXITER
    disp('Did not converge\n');
    fprintf('\t ferr = %e,\t norm(p)=%d\n',ferr,norm(p));
    err=-1;
    estdx=x;
else
    % Two more iterations to polish the estimate
    df=sum(dE2(x,rhat(idx),n,idx));
    d2f=reshape(sum(d2E2(x,rhat(idx),n,idx),1),2,2);
    df=df(:);
    p=-d2f\df;    
    [x,f]=linmin(x,p,func,{rhat,n,idx},1.0e-8,0);        
    fold=f;
    %x=x-d2f\df; %Old code - no line search
    
    df=sum(dE2(x,rhat(idx),n,idx));
    d2f=reshape(sum(d2E2(x,rhat(idx),n,idx),1),2,2);
    df=df(:);
    p=-d2f\df;    
    [x,f]=linmin(x,p,func,{rhat,n,idx},1.0e-8,0);            
    %x=x-d2f\df;%Old code - no line search
     
    estdx=x;
    ferr=abs(f-fold)/max(abs(f),1);
    if debug 
        fprintf('Converged\n');
        fprintf('\t ferr = %e,\t norm(p)=%d\n',ferr,norm(p));
        fprintf('\t estdx=[%5.3f %5.3f]\n',...
            [-x(1) -x(2)]);
        if ~isempty(refdx)
            err=estdx(:)+refdx(:);
            fprintf(' \t refdx=[%5.3f %5.3f]\n \t err=[%5.3e  %5.3e]\n',...
                [refdx(1) refdx(2)],err);
        end
    end
end

