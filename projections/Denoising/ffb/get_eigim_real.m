function [ eigim_real, topeig, topi, topk] = get_eigim_real(U, egval, L0, R, c, n_pc)
% Description
% This function gives mean and denoised images in real space
%
% Input: 
%	U: eigenvectors of C^{(k)}, saved in cell array
% 	fn: inverse Fourier transform of the Fourier- Bessel function in real space on a Cartesian grid.
%	mean_coeff: mean Fourier-Bessel expansion coefficients for k=0
%	Coeff: sPCA expansion coefficients (Wienter filtered).
%	L0: original image size
%	n_max: number of images you want to denoise
%	 
% Output:
%	mean_image: mean image in real space on a Cartesian grid of size L0 x L0
%	denoised: denoised images in real space on a Cartesian grid L0 x L0
% Update 10/15 Zhizhen Zhao

[ fn ] = IFT_FB(R, c);
max_ang_freqs = size(U, 1)-1; %Should be the maximum angular frequency
L = 2*R;
%Computes eigen images, need output from IFT_FB.m.
eig_im = cell(max_ang_freqs+1, 1);
for k = 1:max_ang_freqs+1
    if size(U{k},2)~=0 
        tmp = fn{k};
	tmp = reshape(tmp, L^2, size(tmp, 3));
        tmp2 = tmp*U{k};
        eig_im{k} = reshape(tmp2, L, L, size(tmp2, 2));
    end;
end;

%Original image size
origin = floor(L0/2) + 1;

topeig=[];
topk=[];
topi=[];

for i=1:length(egval)
	bool_pc{i}=zeros(size(egval{i}));
end

for nn=1:n_pc
largest=0;
for i=1:max_ang_freqs+1
	for j=1:length(egval{i})	
		if (egval{i}(j)>largest && bool_pc{i}(j)==0)
			largest=egval{i}(j);
			kval=i; ival=j;
		end	
	end
end
bool_pc{kval}(ival)=1;
topeig=horzcat(topeig,largest);
topk=horzcat(topk,kval);
topi=horzcat(topi,ival);
end

for jj=1:n_pc
	test = zeros(L0);
	test(origin-R:origin+R-1, origin-R:origin+R-1) = real(eig_im{topk(jj)}(:,:,topi(jj)));
	eigim_real(:, :, jj) = test;
end

