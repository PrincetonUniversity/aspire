function [ est_mean ] = mean_LS_batch( ctf, index, regu, nim , nbatch, L, fpath)
%Estimate the mean of the dataset using least squares
% Tejal Bhamre, Oct 2015

ndef=size(ctf,3);
ctf=reshape(ctf,L^2,ndef);
t1=zeros(L^2,1);
t2=zeros(L^2,1);
for nb=1:nbatch
	filename=fullfile(fpath, sprintf('set%d',nb));
	load (filename);
	offset=(nb-1)*(nim/nbatch);
	y=reshape(curr_batch,L^2,size(curr_batch,3));
	
	for i=1:ndef
	    t1=t1+numel(find(index==(i+offset)))*ctf(:,i).^2;
	end
	
	low=(nb-1)*(nim/nbatch)+1;
	high=min(nim,low+(nim/nbatch)-1);
	for i=low:high
	    t2=t2+ctf(:,index(i)).*y(:,i-offset);
	end
end

t1=t1+regu*ones(L^2,1);
est_mean=t2./t1;

est_mean=reshape(est_mean,L,L);
end

