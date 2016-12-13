function mapFile = cryo_fetch_emdID(emdID,verbose)
% CRYO_FETCH_EMDID  Fetch a density map from EMDB
%
% mapFile = cryo_fetch_emdID(emdID)
%   Fetch the map file (MRC format) with the given emdID (integer) from
%   EMDB. The file is downloaded and unzipped.
%   Returns the path to the retrived map file.
%
% mapFile = cryo_fetch_emdID(emdID,verbose)
%   Set verbose to nonzero to print progress messages (default is 1.)
%
% Based on the script
% http://www.zmbh.uni-heidelberg.de/Central_Services/Imaging_Facility/Matlab/Codes/ex6/ex6_8.m
%
% Yoel Shkolnisky, June 2016.

if ~exist('verbose','var')
    verbose=1;
end

currentsilentmode=log_silent(verbose==0);

emdIDstr=int2str(emdID);

ftpServer = 'ftp.ebi.ac.uk';
ftpAddress = ['/pub/databases/emdb/structures/EMD-' emdIDstr '/map/emd_' emdIDstr '.map.gz'];

log_message('Establishing an FTP connection with the EMD server ...');
ngdc = ftp(ftpServer);
log_message('FTP connection was established.');

log_message('Downloading the zipped density map ...');
target=tempname;    % Create a temporary location to save the retrieved map.
zippedMapFile =mget(ngdc, ftpAddress,target);
log_message('The zipped density map was downloaded to %s.',target);

log_message('Closing the FTP connection ...');
close(ngdc);
log_message('The FTP connection was closed.');

log_message('Unzipping the downloaded file ...');
mapFile = gunzip(zippedMapFile{1});
log_message('Unzipping was completed.');
mapFile = mapFile{1};

log_message('(Unzipped) Map file on local computer is : %s',mapFile);

log_silent(currentsilentmode);