symmetry_degree = 8;
rotations_num = 20;
rotations = quat2rotm(randrot(rotations_num, 1));
RelativeRotationEstimationDn(rotations, symmetry_degree);
