function [ est_mean ] = mean_LS( ctf, index, y, regu )
%Estimate the mean of the dataset using least squares
% Tejal, Oct 2015

L=size(y,1);
nim=size(y,3);
ndef=size(ctf,3);

y=reshape(y,L^2,nim);
ctf=reshape(ctf,L^2,ndef);

t1=zeros(L^2,1);
for i=1:ndef
    t1=t1+numel(find(index==i))*ctf(:,i).^2;
end
t1=t1+regu*ones(L^2,1);

t2=zeros(L^2,1);
for i=1:nim
    t2=t2+ctf(:,index(i)).*y(:,i);
end

est_mean=t2./t1;

est_mean=reshape(est_mean,L,L);
end

