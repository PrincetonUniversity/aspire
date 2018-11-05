% INSTALL Install the ASPIRE package
%
% Usage
%    install;
%
% Description
%    Install any external packages and compile all ASPIRE MEX files. This
%    function should only be called once when ASPIRE is first installed.

function install()

    % Set location for temporary files
    tempdir=fmtinput('Enter folder to temporary files ','/tmp','%s');

    try
        fid=fopen(fullfile(aspire_root(),'tmpdir.cfg'),'w');
        fprintf(fid,tempdir);
        fclose(fid);
    catch E
        error('Failed to create tmpdir.cfg Error: %s', E.message);
    end

    try
        tempmrcdir;
    catch E
        error('Failed to create temproray folder. Error: %s',E.message);
    end

    % Install NUFFT
    nufft_choice = multichoice_question( ...
        'Choose a NUFFT package to install ', ...
        {'cims', 'chemnitz', 'finufft', 'none'}, [2 3 4 1], 'cims');

    if nufft_choice == 2
        install_cims_nufft;
    elseif nufft_choice == 3
        install_chemnitz_nfft;
    elseif nufft_choice == 4
        install_finufft;
    end

    % Install SDPLR
    sdplr_choice = two_options_question( ...
        'Install SDPLR? ', 'y', 'n', 'y', '%c');

    if sdplr_choice == 1
        install_sdplr;
    end

    install_mex;
end
