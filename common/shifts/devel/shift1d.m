% Phase correlation shift registration of 1D signals.
%
% Yoel Shkolnisky, December 2013.

B = 1;   % The sampling interval is [-B,B].
N = B*32; % Number of samples in the sampling interval
xx = (-N:1:N)*B/N; %Sampling points in the time domain. Must have odd length.
                   % Pixel size is B/N.
sig2 = (0.1)^2;  % Standard deviation of the sampled Gaussian.
x0 = 0.0; % Center of the Guassisn.
g = @(x) (1/sqrt(2*pi*sig2))*exp(-(x-x0).^2/(2*sig2)); % The Gaussian in the time domain
fg = @(omega) exp(-omega.^2*sig2/2 -sqrt(-1)*x0*omega ); % It's Fourier transform.

s1=g(xx);  % Sample the Guassian
dx=0.2514*B/N;  % We will shift the Guassian by 2 pixels. 
s2=g(xx-dx); % shift to the right.
figure; % Plot original and shifted signals.
plot(xx, s1,'.-');
hold on;
plot(xx,s2,'r.-'); 
hold off;

% Phase correlation
hats1 = fftshift( fft( ifftshift( s1) ))*(B/N); % Compute centered Fourier transform
hats2 = fftshift( fft( ifftshift( s2) ))*(B/N);
rhat=hats1.*conj(hats2)./(abs(hats1.*conj(hats2)));

if any(isnan(rhat)) || any(isinf(rhat))
    error('phase factors have nans or infs');
end

omega = (-N:1:N)*2*pi/B*(N/(2*N+1));
hats1ref=fg(omega);
hats2ref=fg(omega).*exp(-1i.*omega.*dx);

% The following errors should be tiny
fprintf('Err in hats1 = %e\n',norm(hats1-hats1ref));
fprintf('Err in hats2 = %e\n',norm(hats2-hats2ref));

% Plot true and estimated phase factors
figure;
subplot(1,2,1);
plot(real(rhat),'x-b')
hold on;
plot(real(exp(1i.*omega.*dx)),'.-r')
hold off;
subplot(1,2,2);
plot(imag(rhat),'x-b')
hold on;
plot(imag(exp(1i.*omega.*dx)),'.-r')
hold off;

% Verify the expression for the phase factors in terms of N and the number
% of shifted pixels (eliminating B fro the expression)
deltax=dx*N/B; % dx is the shifts in the units of the time domain. 
               % delta x is the number of shifted pixels.
rhat2=exp(2*pi*1i*(-N:N)*deltax/(2*N+1)); % These are the shift pahse factors independent of B.
norm(rhat2-exp(1i.*omega.*dx))

% Estimate shift
r=fftshift(ifft(ifftshift(rhat)));
[vv,ii]=max(r);
c=(numel(xx)-1)/2+1;
est_shift=c-ii;
fprintf('estimated shift = %d\n',est_shift);

% Plot true and estimated phase correlation coefficients
figure;
plot(r,'bx-')
hold on;
rtrue=fftshift(ifft(ifftshift(exp(1i.*omega.*dx))));
plot(rtrue,'r.-')
hold off;


% Optimize for deltax 
%E = @(j,deltax,rho,N) (abs(exp(2*pi*1i*j*deltax/(2*N+1))-rho));
%norm(E(-N:N,deltax,exp(1i.*omega.*dx),N)) % % Test E: result should be tiny

%dE = @(j,deltax,rho,N) ((exp(2.*pi.*1i.*j.*deltax./(2*N+1)).*conj(E(j,deltax,rho,N))-...
 %                      exp(2.*pi.*1i.*j*deltax/(2*N+1)).*E(j,deltax,rho,N))...
 %                      .*2.*pi.*j/(2*N+1));

% Run newton iterations to minimize E. 
% Do not run newton iterations to find the zero of E, as it is not
% applicable in higher dimensions.
MAXITER=10;
eps=1.0e-8;
iter=1;
x=est_shift;
j=-N:N; 
rho=exp(1i.*omega.*dx);

f=sum(abs(E(j,x,rho,N)).^2);
%f=E(x);
while iter<=MAXITER && abs(f)>eps
    df=sum(dE(j,x,rho,N));
    d2f=sum(d2E(j,x,rho,N));
    %d2f1=(sum(dE(j,x+1.0e-4,rho,N))-sum(dE(j,x-1.0e-4,rho,N)))./(2*1.0e-4)
    
    x=x-df/abs(d2f); 
    f=sum(abs(E(j,x,rho,N)).^2);
    iter=iter+1;
end

if iter>=MAXITER
    disp('Did not converge');
else
    % Two more iterations to polish the estimate
    df=sum(dE(j,x,rho,N));
    d2f=sum(d2E(j,x,rho,N));
    x=x-df/abs(d2f);
    df=sum(dE(j,x,rho,N));
    d2f=sum(d2E(j,x,rho,N));
    x=x-df/abs(d2f);
    
    fprintf('Converged on iteration %d, x=%e, err=%e\n',iter-1,x,x-deltax);
end
