function angle_ind = find_angle_mod(angles, value, mod_base, tol)
    angle_ind = find(are_almost_equal(angles, mod(value, mod_base), tol),...
        1, 'first');
    if angle_ind > mod_base + 1e-8
        angle_ind = ceil(angle_ind - mod_base);
    end
end

function result = are_almost_equal(a, b, err)
    result = abs(a - b) <= err;
end
