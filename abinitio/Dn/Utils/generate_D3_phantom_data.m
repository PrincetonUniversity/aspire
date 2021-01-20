function [rots, pf] = generate_D3_phantom_data(n, snr, noise_type, seed)
    vol = make_phantom_D3();
    view3d(vol);
    nr = 89;
    vol = cryo_downsample(vol,nr); % downsample volume
    rots = rand_rots(n); % Generate random rotations
    rots = permute(rots,[2,1,3]); %transpose all matrices
    projs = cryo_project(vol,rots);%project molecule in these orientations
    projs = permute(projs,[2,1,3]);
    viewstack(projs,4,4); % view projections
    rots = permute(rots,[2,1,3]);
    masked_projs = mask_fuzzy(projs,0.5*(nr-1)); % remove noisy edges
    if snr ~= 0 && snr ~= Inf
        masked_projs = cryo_addnoise(masked_projs, snr, noise_type, seed);
    end
    ntheta = 360;
    [pf,~]=cryo_pft(masked_projs,nr,ntheta); % polar fourier transform
end