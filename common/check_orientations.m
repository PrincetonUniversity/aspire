function [h,err_in_degrees,stats]=check_orientations(PHI,PHIref,noplot,verbose)
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
% Set noplot to nonzero to avoid showing the histogram of errors (defult is
% 0, that is plot by default). Set verbose to 0 to eliminate printouts
% (default is 1).
%
% Returns the handle h to the histogram figure, the list of errors between
% each two corresponding Fourier rays, and a struct with statistics.
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
% April 4, 2016 Y.S.    Add noplot and verbose flags.
% Revision:
% Y.S.  November 16,2011. Return the list of estimation errors
% (err_in_degrees).
%  

if ~exist('noplot','var')
    noplot=0;
end

if~exist('verbose','var')
    verbose=1;
end

PHI=register_orientations(PHI,PHIref);
cos_angle = PHI(:,1).*PHIref(:,1) + PHI(:,2).*PHIref(:,2) + PHI(:,3).*PHIref(:,3);
err=acos(cos_angle);
err_in_degrees=err*180/pi;

h=-1;
if ~noplot
    hist(err_in_degrees,50);
    h=gcf;
    set(gca,'FontSize',14);
    %title('Estimation error (in degrees)','FontSize',16);
end

stats.mean_err=mean(err_in_degrees);
stats.med_err=median(err_in_degrees);
stats.std_err=std(err_in_degrees);
stats.min_err=min(err_in_degrees);
stats.max_err=max(err_in_degrees);

if verbose
    log_message('Mean error=%8.6e (degrees)\n',stats.mean_err);
    log_message('std=%8.6e (degrees)\n',stats.std_err);
    log_message('Median error=%8.6e (degrees)\n',stats.med_err);
    log_message('min error=%8.6e (degrees)\n',stats.min_err);
    log_message('max error=%8.6e (degrees)\n',stats.max_err);
end
