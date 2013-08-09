  function offset_zxy = reg_offset_xyz_to_zxy(offset_xyz, dim_xyz)
%|function offset_zxy = reg_offset_xyz_to_zxy(offset_xyz, dim_xyz)
%|
%| convert usual offsets for xyz image ordering
%| to offsets appropriate for zxy image ordering
%|
%| Copyright 2009-5-3, Jeff Fessler, University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

if numel(dim_xyz) ~= 3
	fail('expected dim_zxy to be [3]')
end

if ischar(offset_xyz)
	offset_xyz = penalty_offsets(offset_xyz, dim_xyz);
end

dd_xyz = penalty_displace(offset_xyz, dim_xyz);
dd_zxy = dd_xyz(:, [3 1 2]);
dim_zxy = col(dim_xyz([3 1 2]));
offset_zxy = dd_zxy * [1; cumprod(dim_zxy(1:end-1))];
