function [ aligned_images ] = cryo_align_projections( projs, angle, class_VDM, refl, shifts, k, list_recon, average,global_flip)
% Function for aligning images with its k nearest neighbors to generate
% class averages.
%   Input: 
%       data: LxLxP matrix. P projection images of size LxL pixels.
%       angle: Pxl (l>=k) matrix. Rotational alignment angle
%       class_VDM: Pxl matrix. Nearest neighbor list
%       refl: Pxl matrix. 1: no reflection. 2: reflection
%       FBsPCA_data: Fourier-Bessel steerable PCA data with r_max, UU,
%       Coeff, Mean, Freqs
%       k: number of nearest neighbors for class averages
%       max_shifts: maximum number of pixels to check for shift
%       list_recon: indices for images to compute class averages
%   Output:
%       shifts: Pxk matrix. Relative shifts for k nearest neighbors
%       corr: Pxk matrix. Normalized cross correlation of each image with
%       its k nearest neighbors
%       average: LxLxP matrix. Class averages
%       norm_variance: compute the variance of each class averages.
%
% Zhizhen Zhao Feb 2014

% XXX average is sorted

L=size(projs, 1);

N=floor(L/2);
[x, y]=meshgrid(-N:N, -N:N);
r=sqrt(x.^2+y.^2);

aligned_images=zeros(L,L,length(list_recon),k+1);

%generate grid. Precompute phase for shifts
range=-fix(L/2):fix(L/2);
[omega_x,omega_y]=ndgrid(range,range);
omega_x=-2*pi.*omega_x/L; omega_y=-2*pi.*omega_y/L;
omega_x=omega_x(:); omega_y=omega_y(:);

angle=round(-angle);
angle(angle<0) = angle(angle<0)+360;

angle(angle==360)=0;

%Precompute rotation table, precision for rotation is 1 degree
M=cell(359);
for i=1:359
    M{i}=fastrotateprecomp(L, L,i);
end;

% %Go concurrent
% ps=matlabpool('size');
% if ps==0
%     matlabpool open
% end

for j=1:length(list_recon)
    
    projidx=list_recon(j);
    angle_j=angle(projidx, 1:k); %rotation alignment
 
    refl_j=refl(projidx, 1:k);    %reflections
    
    index = class_VDM(projidx, 1:k); % Indices of nearest neighbor images of image projidx
    
    images=projs(:, :, index); % nearest neighbor images
    
    image1=projs(:, :, projidx);
       
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
    
    images = cat(3, image1, images);
    pf_images=zeros(size(images));
    for i=1:k+1
        pf_images(:, :, i)=cfft2(images(:, :, i));
    end;
    dx=real(shifts(projidx,:)); dy=imag(shifts(projidx,:));
    phase= exp(sqrt(-1)*(omega_x*dx+omega_y*dy));
    pf_images=reshape(pf_images, L^2, k+1);
    pf_images_shift=pf_images.*phase;
    
    
    parfor kk=1:k+1
        aligned_images(:,:,j,kk)=icfft2(reshape(pf_images_shift(:,kk),L,L));
    end

    % If reference averages are given, check that averaging the aligned
    % images indeed results in the give average. The two should be equal to
    % machine precision.
    if exist('average','var')
        flag=1-2*global_flip; % If global_flip is 1 then all averages have 
            % been multiplied by (-1), so multiply back to check accuracy.
            av1=mean(squeeze((aligned_images(:,:,j,:))),3);
        err=norm(av1-flag.*average(:,:,j))/norm(average(:,:,j));
        assert(err<1.0e-15);
    end
    
end

%matlabpool close;