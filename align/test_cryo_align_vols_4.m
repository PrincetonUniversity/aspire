% test_cryo_align_vols_4
%
% Test fast alignment algorithm by aligning density maps
% rotate/translated/reflected relative to, for various symmetry groups. 
%
% Yoel Shkolnisky, April 2021.

% First entry is symmetry group, second is EMD code, third is reported
% model resolution in Angstrom.

clear;

test_densities = {...
{'C1',	 2660, 	3.2},...
{'C1',	 6581,	7.0},...
{'C2',	 0667,	6.2},...
{'C3',	 0731,	2.85},...
{'C4',	 0882,	3.3},...
{'C5',	21376,	2.6},...
{'C7',	11516,	2.38},...
{'C8',	21143,	3.63},...
{'C11',  6458,	4.7},...
{'D2',	30913,	1.93},...
{'D3',	20016,	2.77},...
{'D4',	22462,	2.06},...
{'D7',	 9233,	2.1},...
{'D11',	21140,	3.68},...
{'T',	10835,	1.98},...
{'T',	11267,	4.4},...
{'O',	22658,	1.36},...
{'I',	 3528,	2.7},...
{'I',	22854, 	1.56}};

results=cell(numel(test_densities),3);


for testidx=1:numel(test_densities)
    
    %% Generate two density maps.
    
    test_data=test_densities{testidx};
    symmetry=test_data{1};
    emdid=test_data{2};
    
    symgroup=symmetry(1);

    if numel(symmetry)>1
        symorder=str2double(symmetry(2:end));
    else
        symorder=0;
    end
    
    try
        mapfile=cryo_fetch_emdID(emdid);
        vol=ReadMRC(mapfile);
    catch err
        delete(mapfile);
    end
    delete(mapfile);
    
    [R,~,~]=svd(rand(3)); % Generate random rotation
    volRotated=fastrotate3d(vol,R); % Rotate the reference volume by the random rotation
    volRotated=reshift_vol(volRotated,[5 0 0]);
    %
    % Note: it seems that the estimated shift estdx below is equation to (the
    % negative of) R.'*[0 5 0].', where the flip in the position of the 5 is
    % due to the difference convetions of which coordinate is the x axis
    % between the difference functions. Check if the estdx should R.'*[0 5 0].'
    % or R*[0 5 0].'
    
    %% Visualize the two maps.
    %figure(1); clf; view3d(GaussFilt(cryo_downsample(vol,64),0.1),2.0e-4); title('Reference volume');
    %figure(2); clf; view3d(GaussFilt(cryo_downsample(volRotated,64),0.1),2.0e-4); title('Rotated volume');
    
    %% Align
    verbose=1;
    tic;
    [Rest,estdx,vol2aligned]=cryo_align_vols(symgroup,symorder,vol,volRotated,[],R);
    toc
    
    %figure(3); clf; view3d(vol2aligned,2.0e-4); title('Aligned volume');
    
    % Correlation between two original volumes
    c1=corr(vol(:),volRotated(:));
    fprintf('Correlation between two original volumes %7.4f\n',c1);
    
    % Correlation after alignment
    c2=corr(vol(:),vol2aligned(:));
    fprintf('Correlation between original and aligned volume %7.4f\n',c2);
    
    % The original and aligned volume should have high agreement according to
    % the FSC. FSC may degrade due to sub-pixel misalignment errors.
    %plotFSC(vol,vol2aligned,0.5,1);
    
    % Write results into output table
    results{testidx,1}=testidx;
    results{testidx,2}=symmetry;
    results{testidx,3}=c2;
end