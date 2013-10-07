function [ shifts, corr, average, norm_variance ] = align_main( data, angle, class_VDM, refl, k, max_shifts)
% Function for aligning images with its k nearest neighbors to generate
% class averages.
%   Input: 
%       data: LxLxP matrix. P projection images of size LxL pixels.
%       angle: Pxl (l>=k) matrix. Rotational alignment angle
%       class_VDM: Pxl matrix. Nearest neighbor list
%       refl: Pxl matrix. 1: no reflection. 2: reflection
%       k: number of nearest neighbors for class averages
%       max_shifts: maximum number of pixels to check for shift
%   Output:
%       shifts: Pxk matrix. Relative shifts for k nearest neighbors
%       corr: Pxk matrix. Normalized cross correlation of each image with
%       its k nearest neighbors
%       average: LxLxP matrix. Class averages
%       norm_variance: compute the variance of each class averages.
%
% Zhizhen Zhao Aug 2013

P=size(data, 3);
L=size(data, 1);
l=size(class_VDM, 2);

%Check if the number of nearest neighbors is too large
if l<k
    error('myApp:argChk', 'The number of nearest neighbors is too large. \nIt should be smaller or equal to %d', l)
end;

shifts=zeros(P, k);
corr=zeros(P, k);
average=zeros(L, L, P);
norm_variance=zeros(P, 1);

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

parfor j=1:P
    
    angle_j=angle(j, 1:k); %rotation alignment
 
    refl_j=refl(j, 1:k);    %reflections
    
    index = class_VDM(j, 1:k);
    
    images=data(:, :, index); % nearest neighbor images
    
    image1=data(:, :, j);

    
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
    
    pf1=cfft2(image1);
    pf1=pf1(:);
    
    pf_images=zeros(size(images));
    for i=1:k
        pf_images(:, :, i)=cfft2(images(:, :, i));
    end;
    pf_images=reshape(pf_images, L^2, k);
    
    
    pf=bsxfun(@times, phase, pf1);


    C=pf'*pf_images;
    [corr(j, :), id ] = max(C, [], 1);
    
    pf_images_shift=pf_images.*conj(phase(:, id));
    pf2=[pf_images_shift, pf1];
    variance=var(pf2, [], 2);
    norm_variance(j)=norm(variance, 'fro');
    tmp = mean(pf2, 2);
    tmp = reshape(tmp, L, L);
    average(:, :, j) = icfft2(tmp);

    shifts(j, :)=-shifts_list(id, 1) - sqrt(-1)*shifts_list(id, 2);

end