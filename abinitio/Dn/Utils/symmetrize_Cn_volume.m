function symmetric_vol = symmetrize_Cn_volume(vol, symmetry_degree, ...
    normalize_flag)
    temp_vol = vol;
%     if mod(symmetry_degree, 2) == 1
%        temp_vol = fastrotate3z(vol, 45 / symmetry_degree);
%     end
    symmetric_vol = temp_vol;
    for i=2:symmetry_degree
        symmetric_vol = symmetric_vol + fastrotate3z(temp_vol, ...
            360 * (i - 1)/ symmetry_degree);
    end
    
    if normalize_flag
        symmetric_vol = symmetric_vol / symmetry_degree;
    end
end

