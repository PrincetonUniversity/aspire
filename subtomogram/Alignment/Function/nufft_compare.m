%function [E, T] = nufft_compare(nj,n)
% generate random numbers for the nodes, transform them to random
% nonuniform points in the Fourier domain, then transform them back to the
% original nodes.
% n = 50;
% nj = 1000;
n1 = n;
n2 = n;
n3 = n;
f = zeros(n1+1,n2+1,n3+1);
f(1:end-1,1:end-1,1:end-1) = randn(n1,n2,n3)+1i*randn(n1,n2,n3); % random function values


xj = sort((rand(nj,1)*2-1)*pi); % random points [-pi,pi]
yj = sort((rand(nj,1)*2-1)*pi);
zj = sort((rand(nj,1)*2-1)*pi);

%% Greengard node [-n/2, ..., n/2]
% forward
eps = 1e-16;

iflag = -1;

tic
cg1 = nufft3d2(nj,xj,yj,zj,iflag,eps,n1+1,n2+1,n3+1,f);
TFor_greengard=toc;

% tic
% cg2 = dirft3d2(nj,xj,yj,zj,iflag,n1+1,n2+1,n3+1,f);
% toc

% Compute samples direct
k1=-n1/2:n1/2;
k2=-n2/2:n2/2;
k3=-n3/2:n3/2;
[K1,K2,K3]=ndgrid(k1,k2,k3);
k1=K1(:); clear K1;
k2=K2(:); clear K2;
k3=K3(:); clear K3;
cg3=zeros(nj,1);
tic
for j=1:nj
	x1j=xj(j);
	x2j=yj(j);
	x3j=zj(j);
	cg3(j)=sum(f(:).*exp(-1i*(k1*x1j+k2*x2j+k3*x3j)) );
end 
toc

%EForAbs_greengard = max(max(abs(cg3-cg1)))
%EForAbs_greengard2 = max(max(abs(cg2-cg1)))
EForRel_greengard = max(max(abs(cg3-cg1)))/max(max(abs(cg3)));
%EForRel_greengard2 = max(max(abs(cg2-cg1)))/max(max(abs(cg2)))

% backward
iflag = +1;

tic
fg1 = nufft3d1(nj,xj,yj,zj,cg1,iflag,eps,n1+1,n2+1,n3+1);
TBack_greengard = toc;

% tic
% fg2 = dirft3d1(nj,xj,yj,zj,cg1,iflag,n1+1,n2+1,n3+1);
% toc

% Direct computation
fg3=zeros((n1+1)^3,1);
tic
for j=1:(n1+1)^3
	k1j=k1(j);
	k2j=k2(j);
	k3j=k3(j);
	fg3(j)=sum(cg1.*exp(1i*(k1j*xj+k2j*yj+k3j*zj)) );
end %for
toc

% EBackabs_greengard = max(max(abs(fg1(:)-fg3/nj)))
% EBackabs_greengard2 = max(max(abs(fg1-fg2)))
EBackrel_greengard = max(max(abs(fg1(:)-fg3/nj)))/max(max(abs(fg3/nj)));
% EBackrel_greengard2 = max(max(abs(fg1-fg2)))/max(max(abs(fg2)))

%% Keiner NFFT node [-n/2,n/2-1]

% nj=nj; % number of nodes
% n1=n1; % number of Fourier coefficients in first direction
% n2=n2; % number of Fourier coefficients in second direction
% n3=n3; % number of Fourier coefficients in third direction
N=[n1;n2;n3];

freq=[xj yj zj]/(2*pi); %nodes [-1/2,1/2]

% Initialisation
tic
plan=nfft(3,N,nj); % create plan of class type nfft
plan.x=freq; % set nodes in plan
nfft_precompute_psi(plan); % precomputations
TPrecomp_keiner=toc;

% NFFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhat=f(1:end-1,1:end-1,1:end-1); % Fourier coefficients
fhatv=fhat(:);

% Compute samples with NFFT

plan.fhat=fhatv; % set Fourier coefficients
tic
nfft_trafo(plan); % compute nonequispaced Fourier transform
TFor_keiner = toc;
ck1=plan.f; % get samples


% Compute samples direct
k1=-n1/2:n1/2-1;
k2=-n2/2:n2/2-1;
k3=-n3/2:n3/2-1;
[K1,K2,K3]=ndgrid(k1,k2,k3);
k1=K1(:); clear K1;
k2=K2(:); clear K2;
k3=K3(:); clear K3;
ck2=zeros(nj,1);
tic
for j=1:nj
	x1j=freq(j,1);
	x2j=freq(j,2);
	x3j=freq(j,3);
	ck2(j)=sum( fhatv.*exp(-2*pi*1i*(k1*x1j+k2*x2j+k3*x3j)) );
end 
toc

% Compare results
% EForAbs_keiner = max(abs(ck1-ck2))
EForRel_keiner = max(abs(ck1-ck2))/max(abs(ck2));

% Adjoint NFFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation with NFFT
tic
nfft_adjoint(plan);
TBack_Keiner = toc;
fk1=plan.fhat;

% Direct computation
fk2=zeros(n1*n2*n3,1);
tic
for j=1:n1*n2*n3
	k1j=k1(j);
	k2j=k2(j);
	k3j=k3(j);
	fk2(j)=sum( plan.f.*exp(2*pi*1i*(k1j*freq(:,1)+k2j*freq(:,2)+k3j*freq(:,3))) );
end %for
toc

% Compare results
% EBackAbs_Keiner = max(abs(fk1-fk2));
EBackRel_Keiner = max(abs(fk1-fk2))/max(abs(fk2));

%% Yoel NUFFT and Pseudo polar transform??
% omega = [xj yj zj];
% tic
% cy = nufft_t_3d(f,-omega,'double'); % [-n/2, ..., n/2]
% TFor_Yoel = toc;
% 
% % EForAbs_Yoel = max(abs(cy-cg3))
% EForRel_Yoel = max(abs(cy-cg3))/max(abs(cg3));
% 
% % adjoint
% tic;
% fy1=nufft_3d(cy,omega/(2*pi)*n1,'double',n1); % node in [-n/2, ..., n/2-1]
% TBack_yoel=toc;

% tic
% fy2=slow_nufft_3d(cy,omega/(2*pi)*n1,n1); % omega in [-M/2,M/2]
% t1=toc;


% checkerror(n1,fy2,fy1,cy,1.0e-5,t1,TBack_yoel)
% EBackAbs_Yoel = max(abs(fy1(:)-fk2(:)))
% EBackAbs_Yoel2 = max(abs(fy1(:)-fy2(:)))
% EBackRel_Yoel = max(abs(fy1(:)-fk2(:)))/max(abs(fk2(:)));
% EBackRel_Yoel2 = max(abs(fy1(:)-fy2(:)))/max(abs(fy2(:)))

%% Finufft 
% forward
eps = 1e-16;

iflag = -1;

tic
cf1 = finufft3d2(xj,yj,zj,iflag,eps,f);
TFor_finufft=toc;

% tic
% cg2 = dirft3d2(nj,xj,yj,zj,iflag,n1+1,n2+1,n3+1,f);
% toc

% % Compute samples direct
% k1=-n1/2:n1/2;
% k2=-n2/2:n2/2;
% k3=-n3/2:n3/2;
% [K1,K2,K3]=ndgrid(k1,k2,k3);
% k1=K1(:); clear K1;
% k2=K2(:); clear K2;
% k3=K3(:); clear K3;
% cg3=zeros(nj,1);
% tic
% for j=1:nj
% 	x1j=xj(j);
% 	x2j=yj(j);
% 	x3j=zj(j);
% 	cg3(j)=sum(f(:).*exp(-1i*(k1*x1j+k2*x2j+k3*x3j)) );
% end 
% toc

EForRel_finufft = max(max(abs(cg3-cf1)))/max(max(abs(cg3)));

% backward
iflag = +1;

tic
fi1 = finufft3d1(xj,yj,zj,cg1,iflag,eps,n1+1,n2+1,n3+1);
TBack_finufft = toc;

% tic
% fg2 = dirft3d1(nj,xj,yj,zj,cg1,iflag,n1+1,n2+1,n3+1);
% toc

% Direct computation
% fg3=zeros((n1+1)^3,1);
% tic
% for j=1:(n1+1)^3
% 	k1j=k1(j);
% 	k2j=k2(j);
% 	k3j=k3(j);
% 	fg3(j)=sum(cg1.*exp(1i*(k1j*xj+k2j*yj+k3j*zj)) );
% end %for
% toc

% EBackabs_greengard = max(max(abs(fg1(:)-fg3/nj)))
% EBackabs_greengard2 = max(max(abs(fg1-fg2)))
EBackrel_finufft = max(max(abs(fi1(:)-fg3)))/max(max(abs(fg3)));
% EBackrel_greengard2 = max(max(abs(fg1-fg2)))/max(max(abs(fg2)))


%% Fessler NUFFT not the same with other methods, need to shift because 
% started from 0 to n-1 instead of -n/2 to n/2-1
Nd = [n1 n2 n3];
Jd = Nd/2;
Kd = Nd;
n_shift = [n1/2 n2/2 n3/2];
omega = [xj, yj, zj];

tic
st = nufft_init(omega,Nd,Jd,Kd,n_shift);
TPrecomp_fessler = toc;

tic
cf1 = nufft(f(1:end-1,1:end-1,1:end-1),st); % [-n/2, ..., n/2-1]
TFor_Fessler = toc;

tic
cf2 = dtft(f(1:end-1,1:end-1,1:end-1),omega,'n_shift',n_shift);
toc

% EForAbs_Fessler = max(abs(cf1-ck2)) %-4
% EForAbs_Fessler2 = max(abs(cf1-cf2)) %-4
% EForRel_Fessler = max(abs(cf1-ck2))/max(abs(ck2)); %-6
% % EForRel_Fessler2 = max(abs(cf1-cf2))/max(abs(cf2))
% 
% adjoint
tic
ff1 = nufft_adj(cf1, st);
Tback_fessler=toc;

tic
ff2 = dtft_adj(cf1, omega, Nd, n_shift);
toc

EBackAbs_fessler = max(abs(ff1(:)-fk2(:)))
EBackAbs_fessler2 = max(abs(ff1(:)-ff2(:)))
EBackRel_fessler = max(abs(ff1(:)-fk2(:)))/max(abs(fk2(:)));
EBackRel_fessler2 = max(abs(ff1(:)-ff2(:)))/max(abs(ff2(:)))
%%
E = zeros(3,2);
%E(1,1:2) = [EForRel_Fessler EBackRel_fessler];
E(1,1:2) = [EForRel_greengard EBackrel_greengard];
E(2,1:2) = [EForRel_keiner EBackRel_Keiner];
%E(3,1:2) = [EForRel_Yoel EBackRel_Yoel];
E(3,1:2) = [EForRel_finufft EBackrel_finufft];

T = zeros(3,3);
%T(1,1:3) = [TFor_Fessler Tback_fessler TPrecomp_fessler];
T(1,1:2) = [TFor_greengard TBack_greengard];
T(2,1:3) = [TFor_keiner TBack_Keiner TPrecomp_keiner];
%T(3,1:2) = [TFor_Yoel TBack_yoel];
T(3,1:2) = [TFor_finufft TBack_finufft];
%% compare different method
% Eforabs_GK = max(abs(cg1(:)-ck1(:)));
% fg1_2 = reshape(fg1,n1+1,n2+1,n3+1);
% fg1_2 = fg1_2(1:n1,1:n2,1:n3); % cut the extra ms/2 nodes
% Ebackabs_GK = max(abs(fg1_2(:)-fk1/nj))
% Eforabs_KY = max(abs(ck1-cy))
% Ebackabs_KY = max(abs(fk1-fy1(:)))

%end
