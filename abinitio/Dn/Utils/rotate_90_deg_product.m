function rotated_vol = rotate_90_deg_product(volume, times, axis)
rot_angle = 90;
rotated_vol = volume;
for i=1:times
    if all(axis == [0, 0, 1])
       missing_vals = rotated_vol(1, :, :);
    elseif all(axis == [1, 0, 0])
       missing_vals = rotated_vol(:, :, 1);
    end
    
    rotated_vol = imrotate3(rotated_vol, rot_angle, axis, 'linear', 'crop');
    
    if all(axis == [0, 0, 1])
       rotated_vol(:, 1, :) = flip(missing_vals);
    elseif all(axis == [1, 0, 0])
       rotated_vol(1, :, :) = flip(missing_vals', 2);
    end
end
end
