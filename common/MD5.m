function md5HEXstr=MD5(fname)
%
% MD5 Compute MD5 hash of a file
%
% md5HEXstr=MD5(fname)
%   Compute the MD5 hash of a file named fname. Return the has as a string
%   of hex characters.
%
% Yoel Shkolnisky, April 2017

opt.Input='file';
opt.Format='hex';
opt.Method='MD5';
md5HEXstr=DataHash(fname,opt);
