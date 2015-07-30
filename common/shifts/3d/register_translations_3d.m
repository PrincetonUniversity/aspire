function [estdx,err]=register_translations_3d(vol1,vol2,refdx,debug)
%
% REGISTER_TRANSLATIONS_3D  Estimate relative shift between two volumes.
%
% register_translations_3d(vol1,vol2,refdx)
%   Estimate the relative shift between two volumes vol1 and vol2 to subpixel
%   accuracy. The function uses phase correlation to estimate the relative
%   shift to within one pixel accuray, and then refines the estimation
%   using Newton iterations.
%   
%   Input parameters:
%   vol1,vol2 Two volumes to register. Volumes must be odd-sized.
%   refidx    Two-dimensional vector with the true shift between the images,
%             used for debugging purposes. (Optional)
%
%   Output parameters
%   estidx  A two-dimensional vector of how much to shift vol2 to aligh it
%           with vol1. Returns -1 on error.
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

% Take Fourer transform of both volumes and compute the phase correlation
% factors.
hats1 = fftshift( fftn( ifftshift(vol1))); % Compute centered Fourier transform
hats2 = fftshift( fftn( ifftshift(vol2)));
rhat=hats1.*conj(hats2)./(abs(hats1.*conj(hats2)));
%rhat=hats1.*conj(hats2);
rhat(isnan(rhat))=0;
rhat(isinf(rhat))=0;
if any(isnan(rhat(:))) || any(isinf(rhat(:)))
    error('phase factors have nans or infs');
end

n=size(vol1,1);
ll=fix(n/2);
freqrng=-ll:n-ll-1;

if ~isempty(refdx)
    [omega_x,omega_y,omega_z]=ndgrid(freqrng,freqrng,freqrng);
    omegax=exp(-1i.*2*pi.*omega_x.*refdx(1)/n);
    omegay=exp(-1i.*2*pi.*omega_y.*refdx(2)/n);
    omegaz=exp(-1i.*2*pi.*omega_z.*refdx(3)/n);

    true_phases=omegax.*omegay.*omegaz;
    % The minus is because hats2 is conjugated in rhat , which flips the
    % sign of  its phases.    
    
    % If n is even, then discard the frequnecy which has no conjugate
    % symmetry.
    if mod(n,2)==0
        rhat(:,:,1)=1;
        rhat(:,1,:)=1;
        rhat(1,:,:)=1;
        true_phases(:,:,1)=1;
        true_phases(:,1,:)=1;
        true_phases(1,:,:)=1;
    end

    fprintf('L2 error in estimting rhat = %e\n',norm(rhat(:)-true_phases(:))/norm(true_phases(:)));
end


% Compute the relative shift between the images to to within 1 pixel
% accuracy.
% mm is a window function that can be applied to the volumes before
% computing their relative shift. Experimenting with the code shows that
% windowing does not improve accuracy.
%mm=fuzzymask(n,2,floor(n/4),ceil(n/20));
mm=1; % Windowing does not improve accuracy.
r=fftshift(ifftn(ifftshift(rhat.*mm)));
[~,ii]=max(r(:));

%Find the center
cX=fix(n/2)+1; 
cY=fix(n/2)+1;
cZ=fix(n/2)+1; 

[sX,sY,sZ]=ind2sub(size(r),ii);
est_shift=[cX-sX,cY-sY,cZ-sZ];

if debug
    fprintf('estimated shift = (%d,%d,%d)\n',est_shift(1),est_shift(2),est_shift(3));
end


% Refine shift estimation using Newton iterations

MAXITER=200; % Maximal number of Newton iterations
eps=1.0e-8;  % Error to terminate Newton

if isa(vol1,'single') || isa(vol2,'single') 
    eps=1.0e-5;
end

iter=1;      % Iteration number
x=est_shift; % Initialize Newton from the phase correlation estimated shifts.
x=x(:);
p=1; % Newton step size

% Use only phases close the origin.
radius=fix(n/2)*0.5;
[xx,yy,zz]=ndgrid(freqrng,freqrng,freqrng);
idx=find(xx.^2+yy.^2+zz.^2<radius^2);


% rhat2=hats1.*conj(hats2)./(abs(hats1.*conj(hats2)));
% idx2=find(isnan(rhat2) & isinf(rhat2));
% rhat2(idx2)=0;
% [~,idx2]=sort(abs(rhat(:)));
% idx=idx2(end-round(numel(r)*0.05):end);
% rhat=rhat2;

% Note that we never use the top-most and left-most frequnecies of rhat
% (the frequnecies which have no conjugate-symmetric conuterpart) since in
% such a case, for n even, we would need to modify these values in the
% estimated phases to be real. In other words, for n even, these
% frequnecies do not contain the required information, so ignore them.

% A function to evaluate the L2 error of the estimated shifts.
func=@(x,rhat,n,idx) sum(abs(E3(x,rhat(idx),n,idx)).^2);
lstol=1.0e-5; % Line search tolerance.
f=sum(abs(E3(x,rhat(idx),n,idx)).^2);
ferr=abs(f);
failed=0;
while iter<=MAXITER && ferr>eps &&  norm(p)>eps && ~failed
    df=sum(dE3(x,rhat(idx),n,idx));
    d2f=reshape(sum(d2E3(x,rhat(idx),n,idx),1),3,3);
 
    
    df=df(:);
    p=-d2f\df;    
   
    % Line search. Instead of taking a full Newton step from x in the
    % direction p, find the minimum of f starting from x in the direction
    % p. Note that using line search is not the most efficient way (see
    % Numerical Receipes).
    fold=f;   
    xmin=x;
    fmin=f;
    try
        [xmin,fmin]=linmin(x,p,func,{rhat,n,idx},lstol,0);
    catch
        failed=1;
    end
    x=xmin;
    f=fmin;
        
%     xold=x;
%     lambda=1;    
%     x=x+p; % Since we have a double root at the minimum.
%     f=sum(abs(E3(x,rhat(idx),n,idx)).^2);
%     
%     while f>fold && lambda>1.0e-4
%         lambda=lambda/2;
%         x=xold+lambda*p;
%         f=sum(abs(E3(x,rhat(idx),n,idx)).^2);
%     end

    ferr=abs(f-fold)/max(abs(f),1);
    iter=iter+1;
end

if debug
    fprintf('On iteration %d, ',iter-1);
end

if failed
    estdx=-1;
elseif iter>=MAXITER
    disp('Did not converge\n');
    fprintf('\t ferr = %e,\t norm(p)=%d\n',ferr,norm(p));
    err=-1;
    estdx=x;
else
    % Two more iterations to polish the estimate
    
    if ferr>0
        % There is a case which I don't understand in which I get an error
        % in  brentvec (line 36) "ax, bx, cx must be colinear". To avoid
        % that, I do run polishing iterations if ferr=0.
        
        df=sum(dE3(x,rhat(idx),n,idx));
        d2f=reshape(sum(d2E3(x,rhat(idx),n,idx),1),3,3);
        df=df(:);
        p=-d2f\df;
        [x,f]=linmin(x,p,func,{rhat,n,idx},lstol,0);
        fold=f;
        %x=x-d2f\df; %Old code - no line search
        
        df=sum(dE3(x,rhat(idx),n,idx));
        d2f=reshape(sum(d2E3(x,rhat(idx),n,idx),1),3,3);
        df=df(:);
        p=-d2f\df;
        [x,f]=linmin(x,p,func,{rhat,n,idx},lstol,0);
        %x=x-d2f\df; %Old code - no line search
    end
    
    estdx=x;
    ferr=abs(f-fold)/max(abs(f),1);    
    if debug 
        fprintf('Converged\n');
        fprintf('\t ferr = %e,\t norm(p)=%d\n',ferr,norm(p));
        fprintf('\t estdx=[%5.3f %5.3f %5.3f]\n',...
            [-x(1) -x(2) -x(3)]);
        if ~isempty(refdx)
            err=estdx(:)+refdx(:);
            fprintf(' \t refdx=[%5.3f %5.3f %5.3f]\n \t err=[%5.3e  %5.3e %5.3e]\n',...
                [refdx(1) refdx(2) refdx(3)],err);
        end
    end
end

