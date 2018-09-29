% CRYO_RADIAL_CTF_HANDLE Function handle for radial CTF
%
% Usage
%    h_fun = cryo_radial_ctf_handle(pixel_size, ctf_params);
%
% Input:
%    pixel_size: The pixel size in angstroms.
%    ctf_params: A struct containing the fields:
%          - voltage: Voltage in kV.
%          - defocus: Defocus in angstrom.
%          - spherical_aberration: The spherical aberration constant.
%          - amplitude_contrast: The amplitude constrast constant.
%
% Output
%    h_fun: A function handle accepting one input: the radial frequency
%       between 0 and 1/2, where 1/2 is the Nyquist frequency corresponding to
%       the given pixel size.

function h_fun = cryo_radial_ctf_handle(pixel_size, ctf_params)
    lambda = 12.2643247/sqrt(ctf_params.voltage*1e3 + 0.978466*ctf_params.voltage^2);

    dk = 1/pixel_size;

    k2 = -pi/2*2*lambda*ctf_params.defocus;
    k4 = pi/2*1e7*ctf_params.spherical_aberration*lambda^3;

    h_fun = @(k)( ...
        sqrt(1-ctf_params.amplitude_contrast^2) * ...
            sin(k2*(k*dk).^2+k4*(k*dk).^4) - ...
        ctf_params.amplitude_contrast * ...
            cos(k2*(k*dk).^2+k4*(k*dk).^4));
end
