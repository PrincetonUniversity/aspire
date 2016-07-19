function [proj_CTF, CTF, index]=  add_test_ctf(proj)

%Add CTF to a stack of images proj: the CTF is here is a binary filter, 0.5
%inside a radius and 1 outside.
%Tejal updated Oct 2015

dbstop if error

L=size(proj,1);


% Generate projections with CTF
proj_CTF=zeros(size(proj));
index=ones(size(proj,3),1);

h=ones(L,L);
N=floor(L/2);
[x,y]=meshgrid(-N:N);
x=x/N; y=y/N;
r=x.^2 + y.^2;
h(r<0.5)=0.5*h(r<0.5);

for i=1:size(proj,3)
proj_CTF(:,:,i)=(proj(:,:,i)).*h;% generating data set with CTF (in Fourier space)
end

CTF=h;


