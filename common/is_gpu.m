% IS_GPU Determines whether a GPU is present on the system
%
% Usage
%    b = is_gpu();
%
% Output
%    b: True if a GPU is present, false otherwise.

function b = is_gpu()
    b = exist('gpuDeviceCount') && gpuDeviceCount() > 0;
end
