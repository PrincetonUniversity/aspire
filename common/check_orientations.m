function [h,err_in_degrees]=check_orientations(PHI,PHIref)
%
% Compare the computed orientations to reference orientations.
%
% PHI are the estimated orientation. PHIref are the true orientations.
% Each row in PHI and PHIref corresponds to a point on S2, representing the
% direction vector of some Fourier ray.
%
% The function registers PHI to PHIref and plots the histogram of the
% angles between each pair of corresponding rays.
%
% PHI is returned from cryo_reconstruct. PHIref can be computed by
% euler2S2 using the true Euler angles.
%
% Returns the handle h to the histogram figure.
%
% Example:
%   PHIref=euler2S2('./dataset1/g0ta_movements.txt');
%   check_orientations(PHIc,PHIref);
%
% or
%
%  angles=read_angles_file('./dataset2/g0ta_movements.txt');
%  PHIref=euler2S2(angles);
%  check_orientations(PHIc,PHIref);
%
% Yoel Shkolnisky, August 2008.
%
% Revision:
% Y.S.  November 16,2011. Return the list of estimation errors
% (err_in_degrees).
%  

PHI=register_orientations(PHI,PHIref);
cos_angle = PHI(:,1).*PHIref(:,1) + PHI(:,2).*PHIref(:,2) + PHI(:,3).*PHIref(:,3);
err=acos(cos_angle);
err_in_degrees=err*180/pi;
hist(err_in_degrees,50);
h=gcf;
set(gca,'FontSize',14);
%title('Estimation error (in degrees)','FontSize',16);


fprintf('Mean error=%8.6e (degrees)\n',mean(err_in_degrees));
fprintf('std=%8.6e (degrees)\n',std(err_in_degrees));
fprintf('min error=%8.6e (degrees)\n',min(err_in_degrees));
fprintf('max error=%8.6e (degrees)\n',max(err_in_degrees));

