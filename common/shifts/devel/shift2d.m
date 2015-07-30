% Phase correlation shift registration of 2D images.
%
% Yoel Shkolnisky, December 2013.

B = 1;   % The sampling interval is [-B,B]^2.
N = B*32; % Number of samples in the sampling interval
[xx,yy] = meshgrid((-N:1:N)*B/N,(-N:1:N)*B/N); %Sampling points in the spatial domain. Pixel length/width is B/N.
sig2x = (0.1)^2;  % Standard deviation of the sampled Gaussian.
sig2y = (0.1)^2;
x0 = [0.0, 0.0]; % Center of the Guassisn.
g = @(x,y) (1/(sqrt((2*pi)^2*sig2x*sig2y))).*exp(-(x-x0(1)).^2/(2*sig2x)-(y-x0(2)).^2/(2*sig2y)); % The Gaussian in the spatial domain
fg = @(omegax,omegay) exp(-omegax.^2*sig2x/2 -omegay.^2*sig2y/2-sqrt(-1)*(x0(1)*omegax+x0(2)*omegay)); % It's Fourier transform.

s1=g(xx,yy);  % Sample the Guassian
figure(1);
imagesc(s1);
xlabel('x'); ylabel('y');

dx=[1,0.5]*B/N;  % We will shift the Guassian by 2 pixels. 
s2=g(xx-dx(1),yy-dx(2)); % shift to the right.
figure(1); % Plot original and shifted signals.
imagesc(s2);
xlabel('x'); ylabel('y');

% Phase correlation
hats1 = fftshift( fft2( ifftshift( s1) ))*(B/N)^2; % Compute centered Fourier transform
hats2 = fftshift( fft2( ifftshift( s2) ))*(B/N)^2;

if norm(imag(hats1(:)))/norm(hats1(:))>1.0e-14
    error('Imaginary components too large');
end
hats1=real(hats1); % Note that hats2 is not real due to the shift which introduces pahse.

c1=2*pi/B*(N/(2*N+1));
[omegax,omegay] = meshgrid((-N:1:N)*c1,(-N:1:N)*c1);
hats1ref=fg(omegax,omegay);
hats2ref=fg(omegax,omegay).*exp(-1i.*(omegax.*dx(1)+omegay*dx(2)));

% The following errors should be tiny
fprintf('Err in hats1 = %e\n',norm(hats1-hats1ref));
fprintf('Err in hats2 = %e\n',norm(hats2-hats2ref));

rhat=hats1.*conj(hats2)./(abs(hats1.*conj(hats2)));
if any(isnan(rhat(:))) || any(isinf(rhat(:)))
    error('phase factors have nans or infs');
end


% Plot true and estimated phase factors
figure;
subplot(2,2,1);
imagesc(real(rhat))
title('rhat - real part')
subplot(2,2,2);
imagesc(real(exp(1i.*(omegax.*dx(1)+omegay*dx(2)))))
title('true phases - real part');

subplot(2,2,3);
imagesc(imag(rhat))
title('rhat - imaginary part')
subplot(2,2,4);
imagesc(imag(exp(1i.*(omegax.*dx(1)+omegay*dx(2)))))
title('true phases - imaginary part');

% Verify the expression for the phase factors in terms of N and the number
% of shifted pixels (eliminating B fro the expression)
deltax=dx*N/B; % dx is the shifts in the units of the time domain. 
               % delta x is the number of shifted pixels.
[X,Y]=meshgrid(-N:N,-N:N);
rhat2=exp(2*pi*1i*(X.*deltax(1)+Y.*deltax(2))/(2*N+1)); % These are the shift pahse factors independent of B.
norm(rhat2-exp(1i.*(omegax.*dx(1)+omegay*dx(2))))

% Estimate shift
r=fftshift(ifft2(ifftshift(rhat)));
[vv,ii]=max(r(:));
cX=(size(xx,2)-1)/2+1; %Find the center
cY=(size(xx,1)-1)/2+1; %Find the center

[sY,sX]=ind2sub(size(r),ii);
est_shift=[cX-sX,cY-sY];
fprintf('estimated shift = (%d,%d)\n',est_shift(1),est_shift(2));

% Plot true and estimated phase correlation coefficients
rtrue=fftshift(ifft2(ifftshift(exp(1i.*(omegax.*dx(1)+omegay*dx(2))))));
figure;
subplot(1,2,1)
mesh(r);
subplot(1,2,2)
mesh(rtrue)

% Test dE2
rho=exp(1i.*(omegax.*dx(1)+omegay.*dx(2)));
X0=[1.3,1.2]
DX=[1.0e-4 0];
v1=(sum(sum(abs(E2(X0+DX,rho,N)).^2))-sum(sum(abs(E2(X0-DX,rho,N)).^2)))./(2*norm(DX))
dy=dE2(X0,rho,N);
v2=sum(sum(dy(:,:,1)))
abs(v1-v2)/abs(v1)


% Test d2E2
X0=[1.1,1.7]
dx2=1.0e-4;

% Compute d2x
d2x1=(dE2(X0+[dx2 0],rho,N)-dE2(X0-[dx2 0],rho,N))./(2*dx2);
d2x1=sum(sum(d2x1(:,:,1)))
d2x2=d2E2(X0,rho,N);
d2x2=sum(sum(d2x2(:,:,1)))
errxx=(d2x1-d2x2)./d2x2

% Compute dxdy
% The second coordinate of dE2 is dy which we numerically differentiate by
% x.
dydx1=(dE2(X0+[dx2 0],rho,N)-dE2(X0-[dx2 0],rho,N))./(2*dx2); 
dydx1=sum(sum(dydx1(:,:,2)))
dydx2=d2E2(X0,rho,N);
dydx2=sum(sum(dydx2(:,:,2)))
errxy=(dydx1-dydx2)./dydx2

% Compute d2y
d2y1=(dE2(X0+[0 dx2],rho,N)-dE2(X0-[0 dx2],rho,N))./(2*dx2);
d2y1=sum(sum(d2y1(:,:,2)))
d2y2=d2E2(X0,rho,N);
d2y2=sum(sum(d2y2(:,:,4)))
erryy=(d2y1-d2y2)./d2y2


% Run newton iterations to minimize E
MAXITER=10;
eps=1.0e-8;
iter=1;
x=est_shift;
x=x(:);
j=-N:N; 
rho=exp(1i.*(omegax.*dx(1)+omegay.*dx(2)));
p=1;

f=sum(sum(abs(E2(x,rho,N)).^2));
%f=E(x);
while iter<=MAXITER && abs(f)>eps &&  norm(p)>eps
    df=squeeze(sum(sum(dE2(x,rho,N),1),2));
    d2f=reshape(squeeze(sum(sum(d2E2(x,rho,N),1),2)),2,2);
    %d2f1=(sum(dE(j,x+1.0e-4,rho,N))-sum(dE(j,x-1.0e-4,rho,N)))./(2*1.0e-4)
    
    p=-d2f\df;
    fold=f;
    xold=x;
    lambda=1;    
    x=x+p; % Since we have a double root at the minimum.
    f=sum(sum(abs(E2(x,rho,N)).^2));
    
    while f>fold && lambda>1.0e-4
        lambda=lambda/2;
        x=xold+lambda*p;
        f=sum(sum(abs(E2(x,rho,N)).^2));
    end
    iter=iter+1;
end

if iter>=MAXITER
    disp('Did not converge');
else
    % Two more iterations to polish the estimate
    df=squeeze(sum(sum(dE2(x,rho,N),1),2));
    d2f=reshape(squeeze(sum(sum(d2E2(x,rho,N),1),2)),2,2);
    x=x-d2f\df;
    df=squeeze(sum(sum(dE2(x,rho,N),1),2));
    d2f=reshape(squeeze(sum(sum(d2E2(x,rho,N),1),2)),2,2);
    x=x-d2f\df;
    
    fprintf('Converged on iteration %d, x=[%e %e], err=%e\n',iter-1,[x(1) x(2)],norm(x(:)-deltax(:)));
end
