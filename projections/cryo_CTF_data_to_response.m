function CTFs=cryo_CTF_data_to_response(CTFdata,n,idxlist,scaling)
% CRYO_CTF_DATA_TO_RESPONSE Convert CTF parameters to 2D filter reponse.
%
% CTFs=cryo_phaseflip(CTFdata)   
%   Convert the CTF paramters in CTF data into 2D frequnecy responses.
%   CTFdata contains one record (with CTF parameters) for each projection.
%   It is generated, for example, by reading a STAR file (see example
%   below). CTF records with no corresponding projetctions are ignored. The
%   function returns a 2D nxn freqnuecy response for each record in
%   CTFdata.
%
% CTFs=cryo_CTF_data_to_response(CTFdata,n,idxlist) 
%   Only records of CTFdata whose index is in idxlist are processed. If
%   idxlist is empty or negative, then all records are processed.
% 
% CTFs=cryo_phaseflip(CTFdata,n,idxlist,scaling) 
%   Multiply the pixel size in each CTF record by scaling. This allow to
%   generate the CTF the corresponds to up/down sampled projections.
%
% Example:
%   CTFdata=readSTAR(fname);
%   FPprojs=cryo_CTF_data_to_response(CTFdata,128); 
%
% Yoel Shkolnisky, September 2015.

if nargin<2
    error('CTF dimension not specified');
end

Nctfs=numel(CTFdata.data);
if nargin<3 || isempty(idxlist) || (isscalar(idxlist) && idxlist<0)
    idxlist=1:Nctfs;
end;

if nargin<4
    scaling=1; % No scaling
end
   

CTFs=zeros(n,n,numel(idxlist));

printProgressBarHeader;
for k=1:numel(idxlist)
    idx=idxlist(k);
    progressTic(k,numel(idxlist));
    
    [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
        cryo_parse_Relion_CTF_struct(CTFdata.data{idx});
    h=cryo_CTF_Relion(n,voltage,DefocusU,DefocusV,DefocusAngle,...
        Cs,pixA*scaling,A);

    CTFs(:,:,k)=h;
end
