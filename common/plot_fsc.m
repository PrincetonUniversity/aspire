
function [resA, fighandle] = plot_fsc(fsc, pixelsize, cutoff, fighandle)
%PLOT_FSC get Fourier Shell Correlation curve and plot it,
% along with cutoff line, grid, pixels & angstroms axes, etc.

if ~exist('pixelsize','var'); pixelsize=1.34*359/(numel(fsc)*2+1); end
if ~exist('cutoof','var'); cutoff=0.143; end
if ~exist('fighandle','var'); fighandle=figure; end

n=numel(fsc);

plot(1:n,fsc,'-g','LineWidth',2); % Plot FSC
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

end
