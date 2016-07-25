function ravg3d=cryo_radial_powerspect_3d(vol)
volhat=cfftn(vol);
ravg3d=cryo_radial_average3d((abs(volhat)).^2);
orig3d=ceil((size(ravg3d+1)/2)); % Center of the 3D rdially averaged PSD of the volume.
ravg3d=ravg3d(orig3d(1):end,orig3d(2),orig3d(3)); % Central ray of the PSD.
ravg3d=ravg3d./sum(ravg3d); % Normalize ray to mass 1.