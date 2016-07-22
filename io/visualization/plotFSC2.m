function [resAa,fighandle]=plotFSC2(vol1a,vol2a,vol1b,vol2b,cutoff,pixelsize,fighandle)
%PLOTFSC2 Draw Fourier shell correlation curve and estimate resolution.
%
% Use this function to plot the Fourier shell correlations of two pairs on
% volumes on the same axis. This can be used to compare the quality of two
% reconstructions. 

% plotFSC(vol1a,vol2a,vol1b,vol2b)
%   Draw the Fourier shell correlation curve for the pairs of volumes
%   (vol1a,vol2a) and (vol1b, vol2b), and determine the resolution of each
%   pair according the the given cutoff and pixelsize. Each pair of volumes
%   must be aligned. Uses the default 0.143 cutoff value and a dummy pixels
%   size of 1A.
%
% plotFSC(vol1a,vol2a,vol1b,vol2b,cutoff)
%   Use the given cutoff value for determining resolution. Common values
%   for cutoff are 0.143 and 0.5.  The resolution is determined by the
%   first frequnecy where the correlation goes below this value.
%
% plotFSC(vol1a,vol2a,vol1b,vol2b,cutoff,pixelsize)
%   Use the specified pixel size
%
% plotFSC2(vol1a,vol2a,vol1b,vol2b,cutoff,pixelsize,fighandle)
%   Use the figure corresponding to the given fighandle to plot FSC curve. 
%   If not given, a new handle is allocated.
%
% Output parameters:
%   resA    Resolution of the structure in Angstrom according to the given
%           cutoff value.
%   fighandle   Handle to the FSC curve plot.
%
% Example:
%   vol1a=ReadMRC('./vol1a.mrc');
%   vol2a=ReadMRC('./vol2a.mrc');
%   vol1b=ReadMRC('./vol1b.mrc');
%   vol2b=ReadMRC('./vol2b.mrc');
%   cutoff=0.143;
%   resA=plotFSC(vol1a,vol2a,vol1b,vol2b,cutoff,1.5); 
%
% See plotFSC for more details.
%
% Yoel Shkolnisky, June 2016.



if ~exist('pixelsize','var')
    pixelsize=1;
end

if ~exist('cutoff','var')
    cutoff=0.143;
end

sz1a=size(vol1a);
sz2a=size(vol2a);
sz1b=size(vol1b);
sz2b=size(vol2b);

if any(sz1a-sz2a) || any(sz1a-sz1b) || any(sz1a-sz2b)
    error('All volumes must have te same dimensions');
end

fsc_a=FSCorr(vol1a,vol2a);
n=numel(fsc_a);
fsc_b=FSCorr(vol1b,vol2b);


if ~exist('fighandle','var')
    fighandle=figure;
end


plot(1:n,fsc_a,'-g','LineWidth',2); % Plot FSC

hold on;
plot(1:n,fsc_b,'-r','LineWidth',2); % Plot FSC
hold off;

xlim([1 n]);
ylim([0 1.05]);
grid on

% Plot cutoff line
y=ones(1,n)*cutoff;
hold on; 
plot(1:n,y,'b--','Linewidth',1.5);
hold off;

% Compute resolution - fscres return the bin number where the cutoff
% resolution is obtained.
ja=fscres(fsc_a,cutoff);
resAa=2*pixelsize*n/ja; % Resolution in Angstrom.

jb=fscres(fsc_b,cutoff);
resAb=2*pixelsize*n/jb; % Resolution in Angstrom.


yy=get(gca,'YLim');
line([ja,ja],[0,yy(2)],'LineStyle','--','Color','b','LineWidth',1.5)
line([jb,jb],[0,yy(2)],'LineStyle','--','Color','b','LineWidth',1.5)

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

% Add top axis
assignin('base', 'df', df);
addTopXAxis('XLabStr','Angstroms','expression','1./(argu.*df)');

legend(sprintf('%5.2fA',resAa),sprintf('%5.2fA',resAb));