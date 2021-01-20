run('/scratch/yoelsh/aspire/initpath.m');
symmetry_degree = 4; % n of Dn.
tol = 1e-7;

map = cryo_fetch_emdID(2660); % Download molecule from emd
vol = ReadMRC(map); % Extract volume from .mrc file to 3d array
% At this point use vol to generate D4 volume, and validate that it's D4.
% a = 89;
% vol = reshape(1:a^3, a, a, a);
nr = 89;
vol = cryo_downsample(vol,nr); % downsample volume
tested_vol = symmetrize_D4_volume(vol, true);

for z_angle_index = 1:symmetry_degree - 1
    temp_vol = rotate_90_deg_product(tested_vol, z_angle_index, [0, 0, 1]);
    
    for x_angle_index = [0, 2]
        temp2_vol = rotate_90_deg_product(temp_vol, x_angle_index, [1, 0, 0]);
        assert(all(abs(temp2_vol - tested_vol) < tol, 'all'));
    end
end

log_message('Generation of D4 volume is successfull!');
clear;
