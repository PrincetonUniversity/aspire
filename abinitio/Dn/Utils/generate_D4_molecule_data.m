function [rots, pf] = generate_D4_molecule_data(n, snr, noise_type, seed)
    map = cryo_fetch_emdID(2660); % Download molecule from emd
    vol = ReadMRC(map); % Extract volume from .mrc file to 3d array
    nr = 89;
    vol = cryo_downsample(vol,nr); % downsample volume
    % At this point use vol to generate D4 volume, and validate that it's D4.
%     nr = 3;
%     vol = reshape(1:nr^3, nr, nr, nr);
    vol = symmetrize_D4_volume(vol, true);
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
