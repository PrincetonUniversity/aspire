function masked_projs = circular_mask(projs,c4_type)
%
% Apply a circular mask to each image
% 
% Input parameters:
%   projs      A 3-dimensional array of projection images
%   c4_type    The volume from which the projection images were henerated
%              from. May be 'GAUSSIAN', 'TRPV1' or 'SYNTHETIC'
%
%
% Output parameters:
%   masked_projs  Circulared masked image projections


if strcmp(c4_type,'GAUSSIAN')
    masked_projs = mask_fuzzy(projs,23);
elseif strcmp(c4_type,'TRPV1')
    masked_projs = mask_fuzzy(projs,50);
elseif strcmp(c4_type,'SYNTHETIC')
    masked_projs = mask_fuzzy(projs,50);    
else
    error('implement: find optimal radius for mask');
end

end