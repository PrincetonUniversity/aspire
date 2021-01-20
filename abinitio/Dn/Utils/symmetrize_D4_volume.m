function symmetric_vol = symmetrize_D4_volume(vol, normalize_flag)
    x_rotated_vol = rotate_90_deg_product(vol, 2, [1, 0, 0]);
    symmetric_vol = symmetrize_C4_volume(vol, normalize_flag);
    symmetric_vol = symmetric_vol + ...
        symmetrize_C4_volume(x_rotated_vol, normalize_flag);
    
    if normalize_flag
        symmetric_vol = symmetric_vol / 2;
    end
end
