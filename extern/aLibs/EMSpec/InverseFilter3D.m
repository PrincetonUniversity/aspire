% InverseFilter3D.m
% correct a 3D reconstruction according to the ctfs of the constitutent
% images, using a Wiener filter.

infile='/Users/fred/EMWork/ip3r/IP3R37_trim.mrc';
outfile='/Users/fred/EMWork/ip3r/output.mrc';

defs=[1.2 1.3 1.4 2.0 2.2 2.6];
nims=[500 600 700 800 900 1000];
ndefs=numel(defs);
Bs=200*ones(1,ndefs);  % highly optimistic B-factors

epsi=0.01;  % Wiener filter constant.

lambda=EWavelength(200);

Cs=2;
alpha=.07;

[m s]=ReadMRC(infile);
m=double(m);  % avoid errors in FFT
pixA=double(s.rez)/double(s.nx);
n=size(m,1);

% Add up all the phase-flipped CTFs, in 3D
disp('Making CTFs...');
ct=zeros(n,n,n);
for i=1:ndefs
    ct=ct+nims(i)*abs(CTF3(n,pixA,lambda,defs(i),Cs,Bs(i),alpha));
end;
disp(' done.');
ct=ct/sum(nims);  % Normalize

ct1d=ct(n/2+1:n,n/2+1,n/2+1);

% Construct the Wiener filter function
filt=ct./(ct.^2+epsi);

% Plot it
filt1d=filt(n/2+1:n,n/2+1,n/2+1);
figure(1);
plot((0:n/2-1)'/(n*pixA),[ct1d filt1d]);
legend('Composite CTF','Wiener filter');
xlabel('Spatial frequency, A^{-1}');
drawnow;

% Do the filtering
disp('Filtering...');
fm=real(ifftn(fftn(m).*fftshift(filt)));
disp(' done.');
WriteMRC(fm,pixA,outfile);
figure(2)
ShowSections(fm);
