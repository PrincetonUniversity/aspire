function [rots, pf] = generate_Dn_molecule_data(symmetry_degree, ...
    rots_num, order_2_gen, snr, noise_type, seed, use_yoels)
%GENERATE_MOLECULE_DATA_DN Generation of Dn-symmetric structure and some
%   random 2D projections of this structure. The projections may be
%   generated with some noise.
%
%   Input:
%       symmetry_degree - An integer >= 3.
%       rots_num - An integer, the number of requested projections (N).
%       order_2_gen - Either the matrix g_x=diag([1,-1,-1]) or the matrix
%                     g_y=diag([-1,1,-1]).
%       snr - A floating-point number, which is the signal-to-noise ration.
%             If this number is 0 or Inf, no noise is applied.
%       seed - An integer, for reproducibility.
%       use_yoels - A flag, whether to pick Prof. Shkolnisky's example of a 
%                   D3 structure or not.
%
%   Output:
%       rots - A 3D array of randomly-chosen orientations. The i-th
%              image is given by rots(:,:,i).
%       pf - Polar Fourier-Transformoed projection images. The i-th
%            image is given by pf_images(:,:,i).
%       
%   Written by Elad Eatah May 2021. 

    nr = 89;
    % Select the volume for the simulation.
    if symmetry_degree ~= 3 || ~use_yoels
        %If n is not 3, or n==3 and the requested volume is not Prof.
        % Shkolnisky's D3 volume, a specific volume from EMDB is used.
        
        map = cryo_fetch_emdID(2660); % Download molecule from emd
        vol = ReadMRC(map); % Extract volume from .mrc file to 3d array

        vol = cryo_downsample(vol,nr); % downsample volume
        % Approximate the volume to a Dn-symmetric structure.
        vol = symmetrize_Dn_volume(vol, symmetry_degree, order_2_gen, true);
        tol = 2e-1;  % The artificial Dn object is not so accurate...
    else
        % If n==3 and Prof. Shkolnisky's example is requested, we generate
        % this volume on-demand. It is already a D3 volume.
        vol = make_phantom_D3();
        tol = 1e-14;  % Prof. Shkolnisky's D3 phantom is very accurate...
    end
    view3d(vol);
    
    % Test the proximity of the volume to being a Dn-symmetric structure.
    % This is tested w.r.t both g_x and g_y (these are equivalent iff n is
    % even).
    gx = diag([1,-1,-1]);
    max_err = validate_volume_is_Dn(vol, symmetry_degree, gx, tol);
    fprintf('Dn maximal error w.r.t gx: %f\n', max_err);
    gy = diag([-1,1,-1]);
    max_err = validate_volume_is_Dn(vol, symmetry_degree, gy, tol);
    fprintf('Dn maximal error w.r.t gy: %f\n', max_err);
    
    rots = rand_rots(rots_num); % Generate random rotations
    rots_t = permute(rots,[2,1,3]); %transpose all matrices
    projs = cryo_project(vol,rots_t);%project molecule in these orientations
    projs = permute(projs,[2,1,3]);
    %viewstack(projs,4,4); % view projections
    
    masked_projs = mask_fuzzy(projs,0.5*(nr-1)); % remove noisy edges
    
    % Apply random noise to all images.
    if snr ~= 0 && snr ~= Inf
        masked_projs = cryo_addnoise(masked_projs, snr, noise_type, seed);
    end
    ntheta = 360;
    [pf,~]=cryo_pft(masked_projs,nr,ntheta); % polar fourier transform
end
