function [resolutions,fighandle]=plotFSCmany(vols,cutoff,pixelsize,labels,fighandle)
%PLOTFSCmany Draw Fourier shell correlation curve and estimate resolution.
%
% Use this function to plot the Fourier shell correlations of multiple 
% pairs of volumes on the same axis. 
%
% plotFSCmany(vols)
%   Draw the Fourier shell correlation curve for the pairs of volumes
%   (vols{1},vols{2}), (vols{3},vols{4}),... and determine the resolution
%   of each pair according the the given cutoff and pixelsize. If vols is a
%   contains strings, those are interpreted as names of MRC files.
%   Each pair of volumes must be aligned. Uses the default 0.143 cutoff 
%   value and a dummy pixel size of 1A.
%
% plotFSCmany(vols,cutoff)
%   Use the given cutoff value for determining resolution. Common values
%   for cutoff are 0.143 and 0.5.  The resolution is determined by the
%   first frequnecy where the correlation goes below this value.
%
% plotFSCmany(vols,cutoff,pixelsize)
%   Use the specified pixel size
%
% plotFSCmany(vols,cutoff,pixelsize,labels)
%   Use the given labels in the legend of the FSC plot. Provide labels as
%   cell array of strings
%       plotFSCmany(vols,cutoff,pixA,{'1000','1','1/2','1/4'});
%
% plotFSCmany(vols,cutoff,pixelsize,labels,fighandle)
%   Use the figure corresponding to the given fighandle to plot FSC curve. 
%   If not given, a new handle is allocated.
%
% Output parameters:
%   resoltuions    Resolution of each pair of volumes in Angstrom 
%                  according to the given cutoff value.
%   fighandle   Handle to the FSC curve plot.
%
% Example:
% vols={'volref.mrc','vol_25_1_aligned.mrc',...
%     'volref.mrc','vol_25_2_aligned.mrc',...
%     'volref.mrc','vol_25_3_aligned.mrc',...
%     'volref.mrc','vol_25_4_aligned.mrc'};
% 
% plotFSCmany(vols,cutoff,pixA,{'1000','1','1/2','1/4'});
%
% See plotFSC for more details.
%
% Yoel Shkolnisky, October 2022



if ~exist('pixelsize','var')
    pixelsize=1;
end

if ~exist('cutoff','var')
    cutoff=0.143;
end

if mod(length(vols),2)==1
    error('vols must contain an even number of volumes or filenames')
end

clrs = ['r','g','b','c','m','y','k'];
markers = ['o','+','*','.','x',"square","diamond"];
resolutions = zeros(length(vols)/2,1);
res_labels = cell(length(vols)/2,1);
ja = zeros(length(vols)/2,1);
hplots = zeros(length(vols)/2,1);
sz=-1;

for k=1:length(vols)/2
    
    vol1 = vols{2*k-1};
    vol2 = vols{2*k};
    
    if ischar(vol1) || isstring(vol1)
        vol1=ReadMRC(vol1);
    end
    
    
    if ischar(vol2) || isstring(vol2)
        vol2=ReadMRC(vol2);
    end
    
    if sz==-1
        sz=size(vol1);
    end
    
    sz1=size(vol1);
    sz2=size(vol2);
    
    if any(sz-sz1) || any(sz-sz2)
        error('All input volumes must have te same dimensions');
    end
    
    fsc_a=FSCorr(vol1,vol2);
    n=numel(fsc_a);
    
    if ~exist('fighandle','var')
        fighandle=figure;
    end
    
    if k>1
        hold on;
    end
    hh = plot(1:n,fsc_a,'-'+markers(k)+clrs(k),'LineWidth',2); % Plot FSC
    hplots(k) = hh;

    if k>1
        hold off
    end
       
    % Compute resolution - fscres return the bin number where the cutoff
    % resolution is obtained.
    ja(k)=fscres(fsc_a,cutoff);
    resolutions(k)=2*pixelsize*n/ja(k); % Resolution in Angstrom.
    if ~isempty(labels)
        if length(labels)~=length(vols)/2
            error('Length of labels array must be the same as number of volume pairs');
        end
        str_width = max(strlength(labels));
        formatting_string = strcat('%',int2str(str_width),'s (%5.2fA)');
        res_labels{k}=sprintf(formatting_string,labels{k},resolutions(k));
    else
        res_labels{k}=sprintf('%5.2fA',resolutions(k));
    end
    
end

xlim([1 n]);
ylim([0 1.05]);
grid on

% Plot cutoff line
y=ones(1,n)*cutoff;
hold on;
plot(1:n,y,'b--','Linewidth',1.5);
hold off;
    
yy=get(gca,'YLim');
for k=1:length(vols)/2
    hold on;
    line([ja(k),ja(k)],[0,yy(2)],'LineStyle','--','Color','b','LineWidth',1.5)
    hold off;

end

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

legend(hplots,res_labels);