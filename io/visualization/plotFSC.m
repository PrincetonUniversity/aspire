function [resA,fighandle]=plotFSC(vol1,vol2,cutoff,pixelsize,fighandle)
%PLOTFSC Draw Fourier shell correlation curve and estimate resolution.
%
% plotFSC(vol1,vol2,cutoff,pixelsize)
%   Draw the Fourier shell correlation curve for the volumes vol1 and vol2,
%   and determine the resolution according the the given cutoff and
%   pixelsize. The volumes vol1 and vol2 must be aligned.
%
%   Input parameters:
%   vol1,vol2   3D arrays of the same dimensions. Volumes must be aligned.
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
%   vol1=ReadMRC('./vol1.mrc');
%   vol2=ReadMRC('./vol2.mrc');
%   cutoff=0.143;
%   resA=plotFSC(vol1,vol2,cutoff,1.5); 
%   fprintf('Reolsution from fsc50 curve is %5.2f\n',resA);
%
% Algorithm:
% The resolution resA is determined as follows. 
% The derivation is given in 1D, but it is the same for 3D, if the three
% dimesions of the volume are equal. 
% Denote 
%   B   bandwidth of the signal 
%   fs  Maximal frequnecy 
%   T   Sampling rate (pixel size)
%   N   Length of the signal
%   j   The index of the last frequnecy bin where the correlation is above
%       the cutoff value.
% 
%   1. The bandwidth B is given by B=2*fs, as the signal has frequencies
%      in the range [-fs,fs].
%   2. According to the Nyquist criterion, the samping frequency 1/T must
%      satisfy 1/T>=B. The critical sampling rate is 1/T=B=2*fs.
%   3. The width of each frequnecy bin is B/N=2*fs/N.
%   4. The frequnecy of the bin at index j is fc=2*fs*j/N.
%   5. The resolution is defined as resA=1/fc=N/(2*fs*j).
%      Due to 2 above, resA=N*T/j.
%   6. FSCorr returns correlation values only for positive freqnecuies,
%      that is is returns the correlation at n frequnecies such that N=2*n
%      (it assumes that N is even). Thus, in terms n, resA is given by
%           resA = 2*n*T/j = (2*T)*(n/j).
%
% Yoel Shkolnisky, May 2014.



if ~exist('pixelsize','var')
    pixelsize=1;
end

if ~exist('cutoff','var')
    cutoff=0.143;
end

fsc=FSCorr(vol1,vol2);
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