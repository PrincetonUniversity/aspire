function fighandle=plotPSD(vol,pixelsize,fighandle)
%PLOTPSD Draw Plot power 3D rotational spectrum.
%
% plotPSD(vol,pixelsize)
%   Plot the rotational power spectrum of vol with x-axis units according
%   to the given pixelsize.
%
% plotPSD(vol)
%   Use default pixel size of 1A. 
%
% plotPSD(vol,pixelsize,fighandle)
%   Plot the curve on the given fighandle.
%
% Example:
%   vol=ReadMRC('./vol.mrc');
%   plotPSD(vol,1.34);
%
% Algorithm: See plotFSC
%
% Yoel Shkolnisky, July 2016.



if ~exist('pixelsize','var')
    pixelsize=1;
end

psd=cryo_radial_powerspect_3d(vol);
n=numel(psd);

if ~exist('fighandle','var')
    fighandle=figure;
end

semilogy(1:n,psd,'-g','LineWidth',2); % Plot FSC
xlim([1 n]);
ylim([min(psd)*0.95 1]);
grid on

% Replace the default ticks with frequnecy values
xticks=get(gca,'XTick');
df=1/(2*pixelsize*n);
xticks=xticks*df;

[mjv,mnv]=matlabversion;
if (mjv>8) || (mjv==8 && mnv>=4)
    % Behavior changed on Matlab R2014b
    set(gca,'XTickLabel',sprintf('%7.3f\n',xticks));
else
    set(gca,'XTickLabel',sprintf('%7.3f|',xticks));
end
xlabel('1/A')
ylabel('$\log_{10}(psd)$','Interpreter','Latex');

% Add top axis
assignin('base', 'df', df);
addTopXAxis('XLabStr','Angstrom','expression','1./(argu.*df)');