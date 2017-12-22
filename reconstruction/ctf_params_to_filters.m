% CTF_PARAMS_TO_FILTERS Generate CTF filters from parameters
%
% Usage
%    filter_f = ctf_params_to_filters(L, ctf_params);
%
% Input
%    L: The size of the images.
%    ctf_params: A parameter structure for the CTFs. It is a struct array, with
%       as many elements as CTFs. Each element has the fields
%          - voltage: the voltage of the electron beam in kV,
%          - resolution: angstrom per pixel resolution of the image,
%          - defocus_u: defocus along x-axis in angstrom,
%          - defocus_v: defocus along y-axis in angstrom,
%          - defocus_angle: azimuth of defocus in degrees,
%          - Cs: spherical aberration coefficient, and
%          - alpha: the amplitude contrast.
%
% Output
%    filter_f: An array of size L-by-L-by-K containing the (centered) Fourier
%       transforms of the filters defined by the CTF parameters in ctf_params.

function filter_f = ctf_params_to_filters(L, ctf_params)
    filter_ct = numel(ctf_params);

    filter_f = zeros([L*ones(1, 2) filter_ct]);
    for k = 1:filter_ct
        % NOTE: RELION measures defocus in angstrom while we measure in
        % nanometers.
        filter_f(:,:,k) = -cryo_CTF_Relion(L, ...
            ctf_params(k).voltage, ...
            ctf_params(k).defocus_u/10, ...
            ctf_params(k).defocus_v/10, ...
            ctf_params(k).defocus_angle, ...
            ctf_params(k).Cs, ...
            ctf_params(k).resolution, ...
            ctf_params(k).alpha);
    end
end
