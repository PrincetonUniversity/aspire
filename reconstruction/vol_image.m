% VOL_IMAGE Apply forward imaging model to volume
%
% Usage
%    im = vol_image(vol, params, start, num);
%
% Input
%    vol: A volume of size L-by-L-by-L.
%    params: An imaging parameters structure containing the fields:
%          - rot_matrices: A 3-by-3-by-n array containing the rotation
%             matrices of the various projections.
%          - ctf: An L-by-L-by-K array of CTF images (centered Fourier
%             transforms of the point spread functions of the images)
%             containing the CTFs used to generate the images.
%          - ctf_idx: A vector of length n containing the indices in the
%             `ctf` array corresponding to each projection image.
%          - ampl: A vector of length n specifying the amplitude multiplier
%             of each image.
%    start: The first index of the parameters to use for the imaging mapping.
%    num: The number of images to calculate.
%
% Output
%    im: The images obtained from `vol` by projecting, applying CTFs, and
%       multiplying by the amplitude.

function im = vol_image(vol, params, start, num)
    L = size(vol, 1);

    check_imaging_params(params, L, []);

    range = start:start+num-1;

    im = vol_project(vol, params.rot_matrices(:,:,range));

    im = im_filter(im, params.ctf(:,:,params.ctf_idx(range)));

    im = bsxfun(@times, im, reshape(params.ampl(range), [1 1 num]));

    im = im_translate(im, params.shifts(:,range));
end
