function md5HEXstr=MD5var(data)
%
% MD5 Compute MD5 hash of a variable
%
% md5HEXstr=MD5(data)
%   Compute the MD5 hash of the content of the given variable. 
%   Return the has as a string of hex characters.
%
% Yoel Shkolnisky, April 2017

opt.Input='array';
opt.Format='hex';
opt.Method='MD5';
md5HEXstr=DataHash(data,opt);
