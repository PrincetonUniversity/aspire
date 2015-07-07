function pxx=cryo_noise_estimation_pfs(projs,n_r,n_theta)
% CRYO_NOISE_ESTIMATION_PFS Estimate 2D noise power spectrum in polar
%                                    coordinates.
%
% See Noise_Estimation for more details.
%
% Yoel Shkolnisky, July 2015.


[psd,corr1d]=Noise_Estimation(projs);
R2=cfft2(psd);

ii=norm(imag(R2(:)))/norm(R2(:));
if ii>1.0e-13
    warning('Large imaginary component in R2 = %e',ii);
end
R2=real(R2);

var_n=corr1d(1); % Noise power spectrum
var_s=var(projs(:))-var_n; % Signal power spectrum
snr=var_s/var_n;

psdpf=cryo_pft(R2,n_r,n_theta,'single');
pxx=mean(psdpf,2);
pxx=2*(1/snr)*pxx/mean(pxx); % For compatibility with Mor's code.
