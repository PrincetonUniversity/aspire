function md5HEXstr=MD5(fname)
%
% MD5 Compute MD5 hash of a file
%
% md5HEXstr=MD5(fname)
%   Compute the MD5 hash of a file named fname. Return the has as a string
%   of hex characters.
%
% Yoel Shkolnisky, April 2017

md5HEXstr=GetMD5(fname,'File','hex');
