function [err_in_degrees,mse] = cryo_workflow_Cn_test(xmlname_full)

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

n_symm    = str2double(workflow.algo.n_symm);
n_images  = str2double(workflow.algo.n_images);
snr       = str2double(workflow.algo.snr);
proj_size = str2double(workflow.algo.image_size);
max_shift = ceil(proj_size*str2double(workflow.algo.max_shift_perc)/100);
shift_step = str2double(workflow.algo.shift_step);

[projs,refq,~,~] = generate_cn_images(n_symm,n_images,snr,65,'C1_Eytan',max_shift,shift_step);

[projs,refq] = remove_eq_images(projs,refq);
n_images = size(refq,2);

log_message('loading line indeces cache %s.\n Please be patient...',workflow.cache.name);
load(workflow.cache.name);
log_message('done loading indeces cache');

if snr <= 1
    mask_radius = proj_size*str2double(workflow.algo.mask_radius_perc)/100;
    log_message('Masking images using mask-radius=%d',mask_radius);
    masked_projs = mask_fuzzy(projs,mask_radius);
    
else
    masked_projs = projs;
    log_message('SNR=%.2f is smaller than 1. Not performing mask', snr);
end

precision = 'double';
n_r = ceil(proj_size*str2double(workflow.algo.n_r_perc)/100);
[npf,~] = cryo_pft(masked_projs,n_r,n_theta,precision);

if snr <= 1
    log_message('Guass filtering the images');
    npf = gaussian_filter_imgs(npf);
else
    log_message('SNR=%.2f is smaller than 1. Not performing gauss filtering', snr);
end

log_message('computing self common-line indeces for all candidate viewing directions');
ciis = compute_scls_inds(Ris_tilde,n_symm,n_theta);

is_viz_cls = false;
[vijs,viis,~] = compute_third_row_outer_prod_both_cn(npf,ciis,cijs_inds,Ris_tilde,R_theta_ijs,n_symm,max_shift,shift_step,...
    str2double(workflow.algo.do_handle_equators),refq,is_viz_cls);

vijs = mat2flat(vijs,n_images);
[vijs,viis,~,~] = global_sync_J(vijs,viis);
% 
vis  = estimate_third_rows_ml(vijs,viis);
rots = estimate_inplane_rotations(npf,vis,n_symm,str2double(workflow.algo.inplane_rot_res),max_shift,shift_step);

[rot_alligned,err_in_degrees,mse] = analyze_results_ml(rots,n_symm,n_theta,refq);
%
if str2double(workflow.algo.do_recons_vol)
    log_message('Reconstructing abinitio volume');
    estimatedVol = reconstruct_ml_cn(projs,rot_alligned,n_symm,n_r,n_theta,max_shift,shift_step);
    save_vols(estimatedVol,workflow.algo.recon_mrc_fname,n_symm);
end

close_log();

end


function xmlname_full = create_test_workflow_step1()

workflow_name = '';
while isempty(workflow_name)
    workflow_name = fmtinput('Enter workflow name: ','','%s');
    if ~isvarname(workflow_name)
        fprintf('Workflow name must be a valid filename.\n');
        workflow_name = '';
    end
end

workflow_desc = fmtinput('Enter workflow description: ','','%s');

workflow_dir = fmtinput('Enter full path of output directory: ','','%s');
if ~exist(workflow_dir,'dir') % Do we need to create directory?
    message = 'Output directory does not esxist. Create?';
    do_create = multichoice_question(message,{'Y','N'},[ 1, 0],'Y');
    if do_create == 1
        mkdir(workflow_dir);
    end
else % Check if directory empty
    listing = dir(workflow_dir);
    if numel(listing)>2
        message = 'Output directory not empty. Continue?';
        do_continue = multichoice_question(message,{'Y','N'},[ 1, 0],'N');
        if do_continue == 0
            fprintf('Aborting...\n');
            return;
        end
    end
end

def_logfname = sprintf('log_%s.txt',workflow_name);
workflow_log = fmtinput('Enter log file name: ',def_logfname,'%s');

workflow.info.name = workflow_name;
workflow.info.description = workflow_desc;
workflow.info.created = datestr(now);
workflow.info.working_dir = workflow_dir;
workflow.info.logfile = workflow_log;

open_log(fullfile(workflow.info.working_dir, workflow.info.logfile));

message = 'Rotation symmetry order? (e.g., tpye 3 for C3, type 4 for C4, etc) ';
workflow_n_symm = 1;
while workflow_n_symm < 3
    workflow_n_symm = fmtinput(message,'','%d');
    if workflow_n_symm < 3
        fprintf('rotation symmetry must be an integer equal or greater than 3\n');
    end
end

workflow_n_r_perc = -1;
def_n_r_perc = 50;
while workflow_n_r_perc < 30 || workflow_n_r_perc > 80
    workflow_n_r_perc = fmtinput('Samples in the radial direction (aka n_r) as a percentage of image size? ',def_n_r_perc,'%d');
    if workflow_n_r_perc < 30 || workflow_n_r_perc > 80
       fprintf('Samples in the radial direction percentage should be in the range of [30%%,80%%]\n');
    end
end

do_shifts = multichoice_question('Examine shifts ? ',{'Y','N'},[ 1, 0],'Y');
if do_shifts
    def_per_shift = 15;
    workflow_max_shift_perc = -1;
    while workflow_max_shift_perc <0 || workflow_max_shift_perc > 50
        workflow_max_shift_perc = fmtinput('Maximum shifts to examine in both spatial directions as a percentage of image size? ',def_per_shift,'%d');
        if workflow_max_shift_perc <=0 || workflow_max_shift_perc > 50
            fprintf('Max shift percentage should be in the range of [0%%,50%%]\n');
        end
    end
else
    workflow_max_shift_perc = 0;
end

def_shift_step = 0.5;
if workflow_max_shift_perc > 0
    workflow_shift_step = 0;
    while workflow_shift_step <= 0
        workflow_shift_step = fmtinput('Shift_step? ',def_shift_step,'%d');
        if workflow_shift_step <= 0
            fprintf('Shift_step must be a positive real number\n');
        end
    end 
else
    workflow_shift_step = def_shift_step;
end

clean_snr = 1000000000000;
do_gen_clean_ims = multichoice_question('Generate clean images',{'Y','N'},[ 1, 0],'Y');
if do_gen_clean_ims
    workflow_snr = clean_snr;
else
    workflow_snr = -1;
    while workflow_snr <=0
        workflow_snr = fmtinput('SNR? ',1,'%f');
        if workflow_snr <=0
            fprintf('SNR must be a real positive number\n');
        end
    end
end

if workflow_snr <=1
    do_mask = multichoice_question('Mask images to ignore images corners? ',{'Y','N'},[ 1, 0],'Y');
    if do_mask
        def_per_mask = 40;
        workflow_mask_radius_perc = -1;
        while workflow_mask_radius_perc <=0 || workflow_mask_radius_perc > 80
            workflow_mask_radius_perc = fmtinput('Mask radius as a percentage of image size? ',def_per_mask,'%d');
            if workflow_mask_radius_perc <=0 || workflow_mask_radius_perc > 80
                fprintf('mask radiud should be in the range of [0%%,80%%]\n');
            end
        end
    else
        workflow_mask_radius_perc = 0;
    end
else
    workflow_mask_radius_perc = 0;
end


worflow_do_handle_equators = multichoice_question('handle equator images',{'Y','N'},[ 1, 0],'N');

workflow_inplane_rot_res = 1; % need not be equal to the one used to build the cache
log_message('using inplane rotation resolution of %d degrees', workflow_inplane_rot_res);

workflow_nimages = 0;
while workflow_nimages <=0
    workflow_nimages = fmtinput('Number of images? ',100,'%d');
    if workflow_nimages <=0
        fprintf('Number of images must be a postive number\n');
    end
end

workflow_image_size = 0;
while workflow_image_size <= 0
    workflow_image_size = fmtinput('Size of each image? ',65,'%d');
    if workflow_image_size <= 0 || mod(workflow_image_size,2) == 0
        fprintf('Image size must be a positive odd number');
    end
end

recon_fname = sprintf('%s_c%d_nims%d',workflow_name,workflow_n_symm,workflow_nimages);
recon_mat_fname = strcat(recon_fname,'.mat');
workflow_recon_mat_fname = fullfile(workflow.info.working_dir,recon_mat_fname);
log_message('Reconstructed mat will be written into %s',workflow_recon_mat_fname);
    
worflow_do_recons_vol = multichoice_question('Reconstrcut volume and save to disk? ',{'Y','N'},[ 1, 0],'Y');
if worflow_do_recons_vol
    recon_mrc_fname = strcat(recon_fname,'.mrc');
    workflow_recon_mrc_fname = fullfile(workflow.info.working_dir,recon_mrc_fname);
    log_message('Reconstructed mrc volume will be written into %s',workflow_recon_mrc_fname);
else
    workflow_recon_mrc_fname = ' ';
    log_message('Will not reconstruct volume, nor save it to disk');
end

% workflow.cache.name =  cryo_workflow_abinitio_Cn_cml_reate_cache();

workflow.algo.n_symm             = workflow_n_symm;
workflow.algo.n_images           = workflow_nimages;
workflow.algo.image_size         = workflow_image_size;
workflow.algo.snr                = workflow_snr;
workflow.algo.n_r_perc           = workflow_n_r_perc;
workflow.algo.max_shift_perc     = workflow_max_shift_perc;
workflow.algo.shift_step         = workflow_shift_step;
workflow.algo.mask_radius_perc   = workflow_mask_radius_perc;
workflow.algo.do_handle_equators = worflow_do_handle_equators;
workflow.algo.inplane_rot_res    = workflow_inplane_rot_res;
workflow.algo.do_recons_vol      = worflow_do_recons_vol;
workflow.algo.recon_mrc_fname    = workflow_recon_mrc_fname;
workflow.algo.recon_mat_fname    = workflow_recon_mat_fname;


workflow.cache.name =  cryo_workflow_abinitio_Cn_ml_create_cache();

% Save as xml
tree = struct2xml(workflow);
xmlname = sprintf('%s.xml',workflow_name);
xmlname_full = fullfile(workflow_dir,xmlname);
save(tree,xmlname_full); 

log_message('Workflow file: %s\n',xmlname_full);


close_log();

end