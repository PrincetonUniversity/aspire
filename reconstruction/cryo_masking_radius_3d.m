function [rmin,h1,h2]=cryo_masking_radius_3d(V,cutoff,verbose)
%
% CRYO_MASKING_RADIUS   Estimate masking radius of a projection.
%
% rmin=cryo_masking_radius(P)
%   Automatically estimate the masking radius of the volume V.
%   Find the smallest centered sphere that contains more than 99% of the
%   energy of the volume.
%
% [rmin,h1,h2]=cryo_masking_radius(P,cutoff,verbose)
%   Use the specified cutoff value of the energy. Should be a value between
%   0 and 1. Default is 0.99. 
%   Use cutoff of -1 for default value.
%   Returns handles h1 and h2 to debugging figures.
%
% rmin=cryo_masking_radius(P,cutoff,verbose)
%   Plot progress/debug figures. Default is verbose=0.
%
% Yoel Shkolnisky, October 2016.

if ~exist('verbose','var')
    verbose=0;
end

def_cutoff=0.99;
if ~exist('cutoff','var')
    cutoff=def_cutoff;
end

if cutoff==-1
    cutoff=def_cutoff;
end

if cutoff<0 || cutoff>1
    error('cutoff must be between 0 and 1, or equal to -1 for default value.');
end

if ndims(V)~=3
    error('Input must be a 3D square array');
end

p=size(V,1);
if p~=size(V,2) || p~=size(V,3)
    error('Input must be a 3D square array');
end

c=(p+1)/2; % Center of the image
ravg3d=cryo_radial_average3d(V); % 3D radial average
ravg1d=(squeeze(ravg3d(c+1,c+1,c+1:p))); % 1D radial average
cs=(cumsum(ravg1d.^2)./(norm(ravg1d)^2)); % cs stands for sumsum

if verbose
    log_message('using cutoff=%5.3f',cutoff);
end

rmin=find(cs>cutoff,1,'first');
if verbose
    h1=figure;
    clf;
    plot(cs)
   
    h2=figure;
    clf;
    th = 0:pi/50:2*pi;
    p=size(V,1);
    clr='g';
    
    % X projection
    subplot(1,3,1)
    imagesc(squeeze(sum(V,1)))
    colormap(gray);
    axis image;
    hold on;
    xunit = rmin * cos(th) + c;
    yunit = rmin * sin(th) + c;
    
    % Don't plot parts of the circle that fall outside the image.
    xunit(xunit<=0)=1; xunit(xunit>p)=p;
    yunit(yunit<=0)=1; yunit(yunit>p)=p;
    
    plot(xunit, yunit,'LineWidth',2,'Color',clr);
    hold off;
    title('X')
    
    % Y projection
    subplot(1,3,2)
    imagesc(squeeze(sum(V,2)))
    colormap(gray);
    axis image;
    hold on;
    xunit = rmin * cos(th) + c;
    yunit = rmin * sin(th) + c;
    
    % Don't plot parts of the circle that fall outside the image.
    xunit(xunit<=0)=1; xunit(xunit>p)=p;
    yunit(yunit<=0)=1; yunit(yunit>p)=p;
    
    plot(xunit, yunit,'LineWidth',2,'Color',clr);
    hold off;
    title('Y')
    
    % Y projection
    subplot(1,3,3)
    imagesc(squeeze(sum(V,3)))
    colormap(gray);
    axis image;
    hold on;
    xunit = rmin * cos(th) + c;
    yunit = rmin * sin(th) + c;
    
    % Don't plot parts of the circle that fall outside the image.
    xunit(xunit<=0)=1; xunit(xunit>p)=p;
    yunit(yunit<=0)=1; yunit(yunit>p)=p;
    
    plot(xunit, yunit,'LineWidth',2,'Color',clr);
    hold off;
    title('Z')
end