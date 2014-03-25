load('gaussian_ph.mat');
load p500;

SNR = 1;
[ projections,noisy_projections, defocusID, ctfs ] = ...
    add_ctf_noise( projections,SNR,'gaussian','ctf'  );
fprintf('Reconstruction from clean centered projections affected by CTF\n');
[ v, v_b, kernel ,err, iter, flag] = recon3d_firm_ctf( projections,...
ctfs, defocusID, inv_rot_matrices,[], 1e-6, 1000, zeros(64,64,64));
fprintf('The relative error of reconstruction is %f.\n',...
    norm(v(:)-ph(:))/norm(ph(:)));


[ v, v_b, kernel ,err, iter, flag] = recon3d_firm_ctf_large( projections,...
ctfs, defocusID, inv_rot_matrices,[], 1e-6, 100, zeros(64,64,64),100);
fprintf('The relative error of reconstruction is %f.\n',...
    norm(v(:)-ph(:))/norm(ph(:)));

