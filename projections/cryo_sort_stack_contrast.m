function [idx,vals_sorted]=cryo_sort_stack_contrast(projs)
% CRYO_SORT_STACK_CONTRAST   Sort images by their constrast.
%
% [idx,vals_sorted]=cryo_sort_stack_contrast(projs)
%   Return the indices of the projections sorted by their constrast
%   (variance). 
%
%   Returned values:
%   idx            Indices of sorted images
%   vals_sorted    Constrast (variance) of each image.
%                  vals_sorted(i) is the contrast of image idx(i).
%
% Yoel Shkolnisky, July 2020


averages_contrast=cryo_image_contrast(projs);
[vals_sorted,idx]=sort(averages_contrast,'descend'); 
