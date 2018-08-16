function xmlname_full = cryo_workflow_abinitio_Cn_ml_start

workflow.cache.name =  cryo_workflow_abinitio_Cn_ml_create_cache();

workflow_name='';
while isempty(workflow_name)
    workflow_name = fmtinput('Enter workflow name: ','','%s');
    if ~isvarname(workflow_name)
        fprintf('Workflow name must be a valid filename.\n');
        workflow_name = '';
    end
end

workflow_desc = fmtinput('Enter workflow description: ','','%s');

workflow_mrc_stack_file = '';
while isempty(workflow_mrc_stack_file)
    workflow_mrc_stack_file = fmtinput('Enter full path input MRC stack of images file: ','','%s');
    if exist(workflow_mrc_stack_file,'file')~=2
        fprintf('MRC stack of images file does not exist.\n');
        workflow_mrc_stack_file = '';
    end
end

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
workflow.info.mrc_stack_file = workflow_mrc_stack_file;
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

[~,mrc_header] = ReadMRC(workflow_mrc_stack_file,1, 1);
% Input file has  mrc_header.nz images, each of size mrc_header.nx x
% mrc_header.ny
fprintf('%s :%dx%d (%d projections)\n',...
    workflow_mrc_stack_file, mrc_header.nx, mrc_header.ny, mrc_header.nz);
% message = 'Number of projections to read? ';
% n_images = fmtinput(message,mrc_header.nz,'%d');


workflow_do_downsample = multichoice_question('Downsample images? ',{'Y','N'},[ 1, 0],'N');
if workflow_do_downsample
    workflow_downsample_size = fmtinput('Downsample size in pixels',65,'%d');    
else
    workflow_downsample_size = mrc_header.nx;
end

% do_first_im_ind = multichoice_question('First image ind to select is the first image in stack? ',{'Y','N'},[ 1, 0],'Y');
% if ~do_first_im_ind
%     do_viewstack = multichoice_question('Scroll images before selecting first image ind? ',{'Y','N'},[ 1, 0],'Y');
%     if do_viewstack
%         figure; viewstack(ReadMRC(workflow_mrc_stack_file),5,5);
%     end
%     workflow_first_image_ind = fmtinput('Select first image index',1,'%d');    
% else
%     workflow_first_image_ind = 1;
% end

% workflow_last_image_ind = workflow_first_image_ind + n_images -1;
log_message('images indexes range chosen is [%d,%d]',workflow_first_image_ind,workflow_last_image_ind);

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

do_mask = multichoice_question('Mask images to ignore images corners? ',{'Y','N'},[ 1, 0],'Y');
if do_mask
    def_per_mask = 70;
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
% workflow_mask_radius = fmtinput('Mask_radius in pixels? ',ceil(mrc_header.nx*0.4),'%d');

worflow_do_handle_equators = multichoice_question('handle equator images',{'Y','N'},[ 1, 0],'N');

workflow_inplane_rot_res = 1; % need not be equal to the one used to build the cache
log_message('using inplane rotation resolution of %d degrees', workflow_inplane_rot_res);

recon_fname = sprintf('%s_c%d_ims%dto%d',workflow_name,workflow_n_symm,workflow_first_image_ind,workflow_last_image_ind);
recon_mrc_fname = strcat(recon_fname,'.mrc');
workflow_recon_mrc_fname = fullfile(workflow.info.working_dir,recon_mrc_fname);
log_message('Reconstructed mrc volume will be written into %s',workflow_recon_mrc_fname);

recon_mat_fname = strcat(recon_fname,'.mat');
workflow_recon_mat_fname = fullfile(workflow.info.working_dir,recon_mat_fname);
log_message('Reconstructed mat will be written into %s',workflow_recon_mat_fname);

workflow.algo.n_symm             = workflow_n_symm;
workflow.algo.do_downsample      = workflow_do_downsample;
workflow.algo.downsample_size    = workflow_downsample_size;

% workflow.algo.first_image_ind    = workflow_first_image_ind;
% workflow.algo.last_image_ind     = workflow_last_image_ind;
% workflow.algo.n_images         = workflow_n_images;
workflow.algo.n_r_perc           = workflow_n_r_perc;
workflow.algo.max_shift_perc     = workflow_max_shift_perc;
workflow.algo.shift_step         = workflow_shift_step;
workflow.algo.mask_radius_perc   = workflow_mask_radius_perc;
workflow.algo.do_handle_equators = worflow_do_handle_equators;
workflow.algo.inplane_rot_res    = workflow_inplane_rot_res;
workflow.algo.recon_mrc_fname    = workflow_recon_mrc_fname;
workflow.algo.recon_mat_fname    = workflow_recon_mat_fname;


% Save as xml
tree = struct2xml(workflow);
xmlname = sprintf('%s.xml',workflow_name);
xmlname_full = fullfile(workflow_dir,xmlname);
save(tree,xmlname_full); 

fprintf('Workflow file: %s\n',xmlname_full);
fprintf('Use this file name when calling subsequent funtions.\n');
fprintf('Call next cryo_workflow_abinitio_Cn_ml_execute(''%s'')\n',xmlname_full);

close_log();

end