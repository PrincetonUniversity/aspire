function geo_dist = geodesic_distance(rot1, rot2)
    relative_rot = rot1 * rot2.';
    quat = rotm2quat(relative_rot);
    geo_dist = sqrt(2) * atan2(norm(quat(2:end)), quat(1));
end

