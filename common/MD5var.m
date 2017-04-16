function md5HEXstr=MD5var(data)
%
% MD5 Compute MD5 hash of a variable
%
% md5HEXstr=MD5(data)
%   Compute the MD5 hash of the content of the given variable. 
%   Return the has as a string of hex characters.
%
% Yoel Shkolnisky, April 2017

md5HEXstr=GetMD5(data,'Array','hex');
