function selected=cryo_select_subset(NNtbl,K,toimage,nn_skip)
%
% CRYO_SELECT_SUBSET    Select a subset of the class averages.
%
% selected=cryo_select_subset(NNtbl,K,formimage,toimage)
%   Select a subset of (indices of) K images that are not "too close" to each
%   other. 
%   NNtbl is a nearest neighbors table, where row i contains the
%   indices of the images that are most similar to image i. The rows of the
%   table are assumed to be sort according to "image quality", that is, row
%   1 correspons to the best image (class average), row 2 corresponds to
%   the second best, and so on. The number of rows is the number of images
%   in the dataset.
%   The function selects a subset as follows. It first selects the
%   best-quality image, then discards its neighbors, then selects the best
%   images in the remaining dataset, discards its neighbors, and so on. The
%   number of neighbors discrded at each step is such that we end with a
%   set of K images.
%   If toimage is given, the images are selected only from the range of
%   indices 1...toimage. Since the images are assumed to be sort by
%   quality, this allows choosing images that are "not too bad".
%
%   Input parameters:
%       NNtbl       Nearest neighbors table. A matrix with number of rows
%                   euqal to the number of images in the dataset. Row i in
%                   this table lists the nearest neighbors of image i.
%       K           Number of images to select.
%       toimage     Last image to consider for inclusion in the selected
%                   set. Default: size(NNtbl,1)
%       nn_skip     Number of nearest neighbors to eliminate for each
%                   selected image. Default value is set to allow a final
%                   set of size K.
%
%   Output parameters:
%       selected    List of indices of selected images
%
% Yoel Shkolnisky, August 2018.


N=size(NNtbl,1);    % Number of images in the dataset.

if ~exist('toimage','var')
    toimage=N;
end

if N<K
    error('Cannot select %d images. Only %d images are available.',K,N);
end

if toimage<K
    error('toimage (%d) too small to select %d images',toimage,K);
end

nn_skip1=min(floor(toimage/K),size(NNtbl,2));   % Number of neighbors of 
        % each selected image to eliminate, so that we are rnd up with a
        % subset of K images.
if ~exist('nn_skip','var')
    nn_skip=nn_skip1;
else
    if nn_skip>nn_skip1
        warning('nn_skip too large, setting nn_skip=%d',nn_skip1);
        nn_skip=nn_skip1;
    end
end

mask=ones(toimage,1); % A mask specifying if an image can be selected or 
        % not. If the i element is 0, then image i was a nearest neightbor
        % of a previous image and can not be selected.
        
selected=zeros(K,1); % Indices of selected subset of images.

next2include=1; % Index of the next image to include in the selected subset. 
                % The first image is the one with the highest score.
               
n_selected=1;   % Number of selected images in the subset so far.

while n_selected<=K && next2include<=toimage
    % Select no more than K images, and don't get beyond toimage.
    while ~mask(next2include) && next2include<=toimage
        % Search for the next image that can be selected, but don't get
        % beyond toimage.
        next2include=next2include+1;
    end
    selected(n_selected)=next2include; % The next image to include in the subset.
    
    % Eliminate all neighbors for the currently selected image
    mask(NNtbl(next2include,1:nn_skip))=0;
    next2include=next2include+1;    % The nexy candidate for selection.
    n_selected=n_selected+1;
end

selected=selected(1:n_selected-1);

assert(numel(selected)==K);