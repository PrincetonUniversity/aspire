% Demonstrate the function cryo_estimate_Bfactor

B=zeros(481,1);  % Vector of B estimted factors. 
                 % 481 is the number of micrographs.
fitvals=zeros(size(B));

open_log('Bfactor5.txt');
pixA=1.34; % Pixel size is used to plot axes in units of 1/pixA.
for micrographidx=1:481
    % Load all images of the current micrograph.
    fname=sprintf('%03d_particles_shiny_nb50_new.mrcs',micrographidx');
    inputdir='/mnt/ix2/backup/datasets/80S/pristine/Particles/MRC_1901/';
    projs=ReadMRC(fullfile(inputdir,fname));   
    log_message('************************************************');
    log_message('Loading micrograph #%d (%d projections)',micrographidx,size(projs,3));    
    
    prefix=sprintf('bfactor_%03d',micrographidx);
    [B(micrographidx),fitvals(micrographidx)]=cryo_estimate_Bfactor(projs,pixA,2,prefix,'/home/yoel/tmp/Bfactor');
    
    log_message('B=%7.2f',B(micrographidx));
    log_message('fitvals=%e',fitvals(micrographidx));
end