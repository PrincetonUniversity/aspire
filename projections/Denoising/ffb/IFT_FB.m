% IFT_FB Compute inverse Fourier transform of the Fourier-Bessel basis
%
% Usage
%    fn = IFT_FB(R, c);
%
% Input
%    R: The radius of the disk on which the basis is supported.
%    c: The band limit of the basis in frequency.
%
% Output
%    fn: A cell array containing the inverse Fourier transform (IFT) of the
%       Fourier-Bessel basis. Each cell contains the IFT for a fixed angular
%       frequency, constituted of 2D images arranged along the third dimension,
%       which corresponds to the radial frequency. This is similar to the
%       output of FBcoeff_nfft and recon_images_FB.

% Written by Zhizhen Zhao - 04/2015
% Reformatted, documented, and refactored by Joakim Anden - 2018-Apr-13

function fn = IFT_FB(R, c)
    [x, y] = meshgrid(-R:R-1, -R:R-1);
    r = sqrt(x.^2+y.^2);
    theta = atan2(y, x);

    mask = r(:)<=R;

    theta = theta(mask);
    r = r(mask);

    f = load(fullfile(aspire_root(), 'projections', 'Denoising', 'ffb', ...
        'bessel.mat'));
    bessel = f.bessel(f.bessel(:,4)<=2*pi*c*R,:);
    k_max = max(bessel(:,1));

    fn = cell(k_max+1, 1);
    for i = 1:k_max+1
        bessel_k = bessel(bessel(:,1)==(i-1),:);
        l = size(bessel_k, 1);
        fn_i = zeros((2*R)^2, l);
        for lr = 1:l
            Jk = besselj(i-1, 2*pi*c*r);
            Rkq = bessel_k(lr, 3);
            f_r = Jk./((2*pi*c*r).^2-Rkq^2);
            f_r = 2*c*sqrt(pi)*(-1i)^(i-1)*(-1)^lr*Rkq*f_r;
            f_theta = exp(1i*(i-1)*theta);
            fn_i(mask,lr) = f_theta.*f_r;
        end
        fn{i} = reshape(fn_i, [2*R*ones(1, 2) l]);
    end
end
