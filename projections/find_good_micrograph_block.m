function [r_start,r_end,c_start,c_end,bad]=...
    find_good_micrograph_block(micrograph,min_bad_block,debug)
% FIND_GOOD_MICROGRAPH_BLOCK   Find bad row/columns in a micrograph
%   [...] = find_good_micrograph_block(micrograph,min_bad_block,debug)
%       Find if the micrographs had bad rows/columns. A bad row/column is
%       one whose statistics is significantly different than the statistics
%       of the rest of the micrographs. The function is only looking for
%       bad strips near rhe edges of the micrograph. A bad strip has to be
%       at leat min_bad_block pixels wide to be considered bad. Set debug
%       to non-zero to display the detcted bad area on the micrograph.
%       
% Input parameters:
%   micrographs     A 2D array of the micrograph or a srting with the
%                   filename of the micrograph
%   min_bad_block   A strip of bad pixels must be at least min_bad_block
%                   pixels wide to be considered bad. Default: 10 pixels.
%   debug           Display the micrograph if it contains a bad area, with
%                   the bad area marked. Default: false. 
% 
% Output parameters:
%   r_start     Row where "good" part of the micrograph starts.
%   r_end       Row where "good" part of the micrograph ends.
%   c_start     Column where "good" part of the micrograph starts.
%   c_end       Column where the "good" part of the micrograph ends.
%   bad         Was a bad area detected on the micrograph.
%
% Yoel Shkolnisky, January 2019.


if ~exist('min_bad_block','var') || isempty(min_bad_block)
    min_bad_block = 10; % Minimal number of rows/columns (in pixels) to be 
                        % considered as a damaged area.
end

if ~exist('debug','var')
    debug=false;
end

th=10; % Consider more than th standard deviations as an outlier

if ischar(micrograph)
    micrograph=ReadMRC(micrograph);
end

sz=size(micrograph);

% Find dead areas at the top and bottom
S=sum(micrograph.^2,2)./sz(2);
S_sorted=sort(S);
S_sorted=S_sorted(round(0.2*sz(1)):round(0.8)*sz(1));
mu=mean(S_sorted);
sigma=std(S_sorted);
good=abs(S-mu)<th*sigma;
r_start=find(good,1,'first'); % Start of the "good" zone
r_end=find(good,1,'last'); % End of the "good" zone

% Find dead areas at the left and right
S=sum(micrograph.^2,1)./sz(1);
S_sorted=sort(S);
S_sorted=S_sorted(round(0.2*sz(2)):round(0.8)*sz(2));
mu=mean(S_sorted);
sigma=std(S_sorted);
good=abs(S-mu)<th*sigma;
c_start=find(good,1,'first'); % Start of the "good" zone
c_end=find(good,1,'last'); % End of the "good" zone

bad =  (r_start>min_bad_block) || (r_end<sz(1)-min_bad_block)...
    || (c_start>min_bad_block) || (c_end<sz(2)-min_bad_block);

if debug
    if bad
        p=micrograph;
        % Scale image to be between 0 and 1
        pmin=min(p(:));
        pmax=max(p(:));
        p=(p-pmin)./(pmax-pmin);
        
        im=zeros(sz(1),sz(2),3);
        im(:,:,1)=p;
        im(:,:,2)=p;
        im(:,:,3)=p;
        
        thickness=round(0.01*sz(1));
        if r_start>min_bad_block
            im(r_start:r_start+thickness,:,1)=1;
            im(r_start:r_start+thickness,:,2)=0;
            im(r_start:r_start+thickness,:,3)=0;
        end
        
        if r_end<sz(1)-min_bad_block
            im(r_end-thickness:r_end,:,1)=1;
            im(r_end-thickness:r_end,:,2)=0;
            im(r_end-thickness:r_end,:,3)=0;
        end
        
        if c_start>min_bad_block
            im(:,c_start:c_start+thickness,1)=1;
            im(:,c_start:c_start+thickness,2)=0;
            im(:,c_start:c_start+thickness,3)=0;
        end
        
        if c_end<sz(2)-min_bad_block
            im(:,c_end-thickness:c_end,1)=1;
            im(:,c_end-thickness:c_end,2)=0;
            im(:,c_end-thickness:c_end,3)=0;
        end
        
        imagesc(im);
        axis image;
        drawnow;
    end
end
    