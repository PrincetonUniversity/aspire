function [mjv,mnv,release]=matlabversion
%
% MATALBVERSION     Return MATLAB's version
%
% [mjv,mnv,release]=matlabversion
%   Returns MATLAB's major version, minor version, and release string.
%
% Yoel Shkolnisky, July 2015.

verstr=version; % Get MATLAB's version string
verdata=sscanf(verstr,'%d.%d.%d.%d %s');

mjv=verdata(1); % Major version
mnv=verdata(2); % Minor version
verdata=verdata(5:end); % Strip version codes so we a left with release data.
verdatastr=sprintf('%s',verdata);
releasestartidx=strfind(verdatastr,'('); % Find opening '(';
releaseendidx=strfind(verdatastr,')'); % Find closing '(';
release=verdatastr(releasestartidx(1)+1:releaseendidx(end)-1);
