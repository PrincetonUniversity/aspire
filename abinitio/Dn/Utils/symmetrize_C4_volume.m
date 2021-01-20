function symmetric_vol = symmetrize_C4_volume(vol, normalize_flag)
    symmetry_degree = 4;
    symmetric_vol = vol;
    temp_vol = vol;
    axis = [0, 0, 1];
    for i=1:symmetry_degree - 1
        temp_vol = rotate_90_deg_product(temp_vol, 1, axis);
        symmetric_vol = symmetric_vol + temp_vol;
    end
    
    if normalize_flag
        symmetric_vol = symmetric_vol / symmetry_degree;
    end
end
