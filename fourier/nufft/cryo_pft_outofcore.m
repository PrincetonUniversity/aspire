function freqs=cryo_pft_outofcore(instack,outstack,n_r,n_theta)
%
% Compute the polar Fourier transform of projections with resolution n_r in
% the radial direction and resolution n_theta in the angular direction.
%
% Input parameters:
%   instack    Name of MRC file containing the projections to transform.
%              Images must be real-valued.
%   outstack   Name of MRC file with Fourier transforms of the projections
%              (complex-valued). Use imagestackReaderComplex to read this
%              file.
%   n_r        Number of samples along each ray (in the radial direction).
%   n_theta    Angular resolution. Number of Fourier rays computed for each
%              projection.
%   precision  'single' or 'double'. Default is 'single'. The polar Fourier
%              samples for 'single' are computed to accuracy of 1.e-6.
%
% Output parameters:
%   freqs   Frequencies at which the polar Fourier samples were computed. A
%           matrix with n_rxn_theta rows and two columns (omega_x and
%           omega_y).
%          
% Yoel Shkolnisky, April 2017.

% NOTE: Only sigle precision is supported.
precision='single';

   
%n_uv=size(p,1);
omega0=2*pi/(2*n_r-1);
dtheta=2*pi/n_theta;

freqs=zeros(n_r*n_theta,2); % sampling points in the Fourier domain
for j=1:n_theta
    for k=1:n_r
        freqs((j-1)*n_r+k,:)=[(k-1)*omega0*sin((j-1)*dtheta),...
            (k-1)*omega0*cos((j-1)*dtheta)];
    end
end

%   freqs is the frequencies on [-pi,pi]x[-pi,pi] on which we sample the
%   Fourier transform of the projections. An array of size n_r*n_theta by 2
%   where each row corresponds to a frequnecy at which we sample the
%   Fourier transform of the projections. The first column is omega_x, the
%   second is omega_y. 

imreader=imagestackReader(instack);
% precomputed interpolation weights once for the give polar grid. This is
% used below for computing the polar Fourier transform of all slices
n=imreader.dim(1);
n_projs=imreader.dim(3);

% Get the current parallel pool to query for the number of available
% workers. If no pool exists, create one.
cp=gcp;
nWorkers=cp.NumWorkers;

% Create temporary file to hold the partial PFT computed by each of the
% workers.
fnames=cell(nWorkers,1);
for worker=1:nWorkers
    fnames{worker}=tempmrcname;
end

% Each worker processes chuncksize images.
chuncksize=ceil(n_projs/nWorkers);

parfor worker=1:nWorkers
    % Determine the indices of the images to be processed by the current
    % worker.
    idx=(worker-1)*chuncksize+1:min(chuncksize*worker,n_projs);    
    
    pf=imagestackWriterComplex(fnames{worker},numel(idx),100);
    for k=1:numel(idx)
        tmp=imreader.getImage(idx(k));
        tmp = nufft2(tmp, -freqs');  
        pf.append(reshape(tmp,n_r,n_theta));   
    end
    pf.close;
end

% Merge all temporary files into a single file.
pf=imagestackWriterComplex(outstack,n_projs,100);
for worker=1:nWorkers
    stackreader=imagestackReaderComplex(fnames{worker});
    for k=1:stackreader.dim(3)
        fim=stackreader.getImage(k);
        pf.append(fim);
    end
end
pf.close;

% Delete temporary files
for worker=1:nWorkers
    delete(fnames{worker});
end

