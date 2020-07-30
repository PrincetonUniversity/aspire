function [voltage,DefocusU,DefocusV,DefocusAngle,Cs,pixA,A]=...
    cryo_parse_Relion_CTF_struct(CTFdata,k)
% Tejal Oct 2015, based on Yoel's function
% Parse .star file to extract parameters in RELION format

if numel(CTFdata)==1 % RELION < 3.1
    voltage=CTFdata.data{k}.rlnVoltage;
    DefocusU=CTFdata.data{k}.rlnDefocusU/10; % Relion uses Angstrom. Convert to nm.
    
    if isfield(CTFdata.data{k},'rlnDefocusV')
        DefocusV=CTFdata.data{k}.rlnDefocusV/10; % Relion uses Angstrom. Convert to nm.
    else
        DefocusV=DefocusU;
    end
    
    if isfield(CTFdata.data{k},'rlnDefocusAngle')
        DefocusAngle=CTFdata.data{k}.rlnDefocusAngle*pi/180; % Convert to radians.
    else
        DefocusAngle=0;
    end
    
    Cs=CTFdata.data{k}.rlnSphericalAberration; % In mm, No conversion is needed.
    
    pixA=-1;
    if isfield(CTFdata.data{k},'rlnDetectorPixelSize')
        PS=CTFdata.data{k}.rlnDetectorPixelSize; % In microns. Convert to Angstroms below.
        mag=CTFdata.data{k}.rlnMagnification;
        pixA=PS*10^4/mag; % Convert pixel size on the detector in microns to spatial resoution in Angstroms.
    elseif isfield(CTFdata.data{k},'pixA')
        pixA=CTFdata.data{k}.pixA;
        % else
        %     errmsg=sprintf(strcat('Cannot get pixel size from CTF data.\n',...
        %         'Either pixA or rlnDetectorPixelSize together with rlnMagnification should be given.\n',...
        %         'Use addfieldtoSTARdata to add pixA manually to the STAR file.'));
        %     error('%s',errmsg);
    end
    A=CTFdata.data{k}.rlnAmplitudeContrast;
else
    
    % NOTE: If STAR files of RELION 3.1 is used, then the structure of the
    % STAR file is assumed to contained one optics group (location 1 in the
    % stardata array) and one particles group (location 2 in the stardata
    % array).
    
    voltage=CTFdata(1).data{1}.rlnVoltage;
    DefocusU=CTFdata(2).data{k}.rlnDefocusU/10; % Relion uses Angstrom. Convert to nm.
    
    if isfield(CTFdata(2).data{k},'rlnDefocusV')
        DefocusV=CTFdata(2).data{k}.rlnDefocusV/10; % Relion uses Angstrom. Convert to nm.
    else
        DefocusV=DefocusU;
    end
    
    if isfield(CTFdata(2).data{k},'rlnDefocusAngle')
        DefocusAngle=CTFdata(2).data{k}.rlnDefocusAngle*pi/180; % Convert to radians.
    else
        DefocusAngle=0;
    end
    
    Cs=CTFdata(1).data{1}.rlnSphericalAberration; % In mm, No conversion is needed.
    
    pixA=-1;
    if isfield(CTFdata(1).data{1},'rlnImagePixelSize') % In angstroms
        pixA=CTFdata(1).data{1}.rlnImagePixelSize;
    elseif isfield(CTFdata,'pixA')
        pixA=CTFdata(2).data{k}.pixA;
        % else
        %     errmsg=sprintf(strcat('Cannot get pixel size from CTF data.\n',...
        %         'Either pixA or rlnDetectorPixelSize together with rlnMagnification should be given.\n',...
        %         'Use addfieldtoSTARdata to add pixA manually to the STAR file.'));
        %     error('%s',errmsg);
    end
    A=CTFdata(1).data{1}.rlnAmplitudeContrast;
end