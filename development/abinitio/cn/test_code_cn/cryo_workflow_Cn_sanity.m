function [err_in_degrees,mse] = cryo_workflow_Cn_sanity(xmlname_full)
% cryo_workflow_Cn_sanity('/home/gabip/matlabProjects/aspire/aspire/development/abinitio/C5/ml_cn/test/test_all.xml');
if exist('xmlname_full','var')
    try
        validate_xml_step1(xmlname_full);
    catch
        error('Problem with xml file. Execute this file without any params\n');
    end
else
    xmlname_full = create_test_workflow_step1();
    validate_xml_step1(xmlname_full); % self consistancy check. this should be ok
end

% Read workflow file
tree = xmltree(xmlname_full);
workflow = convert(tree);

open_log(fullfile(workflow.info.working_dir, workflow.info.logfile));

n_symms = 3:4;
err_in_degrees_all = cell(1,numel(n_symms));
mse_all = cell(1,numel(n_symms));
counter = 1;
for n_symm=n_symms
    workflow.algo.n_symm = n_symm;
    [err_in_degrees,mse] = cryo_workflow_Cn_test(xmlname_full);
    err_in_degrees_all{counter} = err_in_degrees;
    mse_all{counter} = mse;
    counter  = counter + 1;
end

close_log();

end