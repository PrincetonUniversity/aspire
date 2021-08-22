function [resA,fighandle]=plotFRC(im1,im2,cutoff,pixelsize,fighandle)
%PLOTFRC Draw Fourier ring correlation curve and estimate resolution.
%
% plotFSC(im1,im2,cutoff,pixelsize)
%   Draw the Fourier ring correlation curve for the images im1 and im2,
%   and determine the resolution according the the given cutoff and
%   pixelsize. The images im1 and im2 must be aligned.
%
%   Input parameters:
%   im1,im2   2D images of the same dimensions. Images must be aligned.
%   cutoff  (Optional) Correlation cutoff threshold to determine
%           resolution. The resolution is determined by the first frequnecy
%           where the correlation goes below this value. Common values are
%           0.143 and 0.5. Default is 0.143.
%   pixelsize (Optional) Default=1A.
%   fighandle Where to plot the figure. If not given, a new handle is
%             allocated.
%
%   Output parameters:
%   resA    Resolution of the structure in Angstrom according to the given
%           cutoff value.
%   fighandle   Handle to the FSC curve plot.
%
% Example:
%   im1=ReadMRC('./vol1.mrc');
%   im2=ReadMRC('./vol2.mrc');
%   cutoff=0.143;
%   resA=plotFSC(im1,im2,cutoff,1.5); 
%   fprintf('Reolsution from fsc curve is %5.2f\n',resA);
%
% Algorithm: See plotFSC.m
% Adapted from plotFSC.m
%
% Yoel Shkolnisky, August 2021.



if ~exist('pixelsize','var')
    pixelsize=1;
end

if ~exist('cutoff','var')
    cutoff=0.143;
end

fsc=FRCorr(im1,im2);
n=numel(fsc);

if ~exist('fighandle','var')
    fighandle=figure;
end

plot(1:n,fsc,'-g','LineWidth',2); % Plot FSC
xlim([1 n]);
ylim([-0.1 1.05]); % To allow presenting oscilations around zero, as suggested by Van Heel.
grid on

% Plot cutoff line
y=ones(1,n)*cutoff;
hold on; 
plot(1:n,y,'b--','Linewidth',1.5);
hold off;

% Compute resolution - fscres return the bin number where the cutoff
% resolution is obtained.
j=fscres(fsc,cutoff);
resA=2*pixelsize*n/j; % Resolution in Angstrom.

yy=get(gca,'YLim');
line([j,j],[0,yy(2)],'LineStyle','--','Color','b','LineWidth',1.5)

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
addTopXAxis('XLabStr',sprintf('Resolution=%5.2fA',resA),'expression','1./(argu.*df)');