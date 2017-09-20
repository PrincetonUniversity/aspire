function [ shifts, corr, averagesfname, norm_variance ] = align_main( data, angle, class_VDM, refl, sPCA_data, k, max_shifts, list_recon, tmpdir)
% Function for aligning images with its k nearest neighbors to generate
% class averages.
%   Input: 
%       data: LxLxP matrix. P projection images of size LxL pixels.
%       angle: Pxl (l>=k) matrix. Rotational alignment angle
%       class_VDM: Pxl matrix. Nearest neighbor list
%       refl: Pxl matrix. 1: no reflection. 2: reflection
%       sPCA_data: steerable PCA data with eig_im, Coeff, Freqs, R
%       Coeff, Mean, Freqs
%       k: number of nearest neighbors for class averages
%       max_shifts: maximum number of pixels to check for shifts
%       tmpdir  temporary folder for intermetidate files. Must be empty.
%   Output:
%       shifts: Pxk matrix. Relative shifts for k nearest neighbors
%       corr: Pxk matrix. Normalized cross correlation of each image with
%       its k nearest neighbors
%       averagesfname: Class averages folder name
%       norm_variance: compute the variance of each class averages.
%
% Zhizhen Zhao Feb 2014
% Tejal Bhamre, 3/2017: Use recon_spca
% Zhizhen Zhao, shouldn't use recon_spca, too much memory used, should compute on the fly
P=data.dim(3);
L=data.dim(1);
l=size(class_VDM, 2);

N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);

r_max = sPCA_data.R;

%Check if the number of nearest neighbors is too large
if l<k
    error('myApp:argChk', 'The number of nearest neighbors is too large. \nIt should be smaller or equal to %d', l)
end;

shifts=zeros(length(list_recon), k+1);
corr=zeros(length(list_recon), k+1);
norm_variance=zeros(length(list_recon), 1);

%generate grid. Precompute phase for shifts
range=-fix(L/2):fix(L/2);
[omega_x,omega_y]=ndgrid(range,range);
omega_x=-2*pi.*omega_x/L; omega_y=-2*pi.*omega_y/L;
omega_x=omega_x(:); omega_y=omega_y(:);
a=-max_shifts:1:max_shifts; % checking shift for every other pixel;
num=length(a);
a1=repmat(a', num, 1);
a2=kron(a', ones(num, 1));
shifts_list=[a1, a2];

phase= exp(sqrt(-1)*(omega_x*shifts_list(:, 1)'+omega_y*shifts_list(:, 2)'));

angle=round(-angle);
angle(angle<0) = angle(angle<0)+360;

angle(angle==360)=0;

%Precompute rotation table, precision for rotation is 1 degree
M=cell(359);
for i=1:359
    M{i}=fastrotateprecomp(L, L,i);
end;


filelist=dir(fullfile(tmpdir,'*.*'));
if numel(filelist)>2
    error('Directory %s is not empty. Aborting',tmpdir)
end

printProgressBarHeader;
parfor j=1:length(list_recon)
    progressTic(j,length(list_recon));
    
    angle_j=angle(list_recon(j), 1:k); %rotation alignment
 
    refl_j=refl(list_recon(j), 1:k);    %reflections
    
    index = class_VDM(list_recon(j), 1:k);
    
    images=data.getImage(index); % nearest neighbor images
    
    image1=data.getImage(list_recon(j));
    
    for i=1:k
        if (refl_j(i)==2)
            images(:, :, i)=flipud(images(:, :, i));
        end;
    end;
    
    for i=1:k
        if (angle_j(i)~=0)
            images(:, :, i)=fastrotate(images(:, :, i), angle_j(i), M{angle_j(i)});
        end
    end;
    
    %Build denoised images from fast steerable PCA
    origin = floor(L/2) + 1;
    R = sPCA_data.R;
    mean_im = sPCA_data.fn0*sPCA_data.Mean;    
    tmp = sPCA_data.eig_im(:, sPCA_data.Freqs == 0)*sPCA_data.Coeff(sPCA_data.Freqs == 0, j) + 2*real(sPCA_data.eig_im(:, sPCA_data.Freqs~=0)*sPCA_data.Coeff(sPCA_data.Freqs~=0, j));
    tmp = tmp + mean_im;
    I = zeros(L);
    I(origin-R:origin+R-1, origin-R:origin+R-1, :) = reshape(tmp, 2*R, 2*R); 
    
    pf1 = cfft2(I);
    pf1 = pf1(:);
    
    images = cat(3, image1, images);
    pf_images=zeros(size(images));
    for i=1:k+1
        pf_images(:, :, i)=cfft2(images(:, :, i));
    end;
    pf_images=reshape(pf_images, L^2, k+1);
    
    
    pf=bsxfun(@times, phase, pf1);

    C=pf'*pf_images;
    [corr(j, :), id ] = max(C, [], 1);
    
    pf_images_shift=pf_images.*conj(phase(:, id));
    variance=var(pf_images_shift, [], 2);
    norm_variance(j)=norm(variance, 'fro');
    tmp = mean(pf_images_shift, 2);
    tmp = reshape(tmp, L, L);
    
    mrcname=sprintf('average%d.mrc',j);
    mrcname=fullfile(tmpdir,mrcname);
    average = real(icfft2(tmp));
    WriteMRC(average,1, mrcname);    
    shifts(j, :)=-shifts_list(id, 1) - sqrt(-1)*shifts_list(id, 2);

end

% Merge all averages into a single file.
averagesfname=tempname;
[~, averagesfname]=fileparts(averagesfname);
averagesfname=fullfile(tmpdir,averagesfname);
stack=imagestackWriter(averagesfname,numel(list_recon),1,100);
for j=1:length(list_recon)
    mrcname=sprintf('average%d.mrc',j);
    mrcname=fullfile(tmpdir,mrcname);
    average = ReadMRC(mrcname);
    stack.append(average);
end
stack.close;

end
