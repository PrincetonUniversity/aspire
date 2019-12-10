function mgraph=micrograph(micrograph_size,particle_size,n_particles)
%
% MICROGRAPH    Simulate clean micrograph
%
% mgraph=micrograph(micrograph_size,particle_size,n_particles)
%   Generate an micrograph of size micrograph_size x micrograph_size
%   containing n_particles, each of size particle_size x particle_size.
%
% Example: mm=micrograph(1024,65,50);
%
% Yoel Shkolnisky, May 2017

if particle_size>micrograph_size
    error('Cannot fit any particle into micrograph');
end

M=floor(micrograph_size/particle_size); % how many particles fit into the 
    % micrograph along each dimension.

if n_particles>M*M
    error('Cannot fit %d particles into micrograph of size %dx%d',...
        n_particles,micrograph_size,micrograph_size);
end
    
mgraph=zeros(micrograph_size);
idx=randperm(M*M); % permute the indices of the blocks where we can put particles.
idx=idx(1:n_particles); % but keep only nparticles indices.

% Generate clean images
dummysnr=1;
max_shift=0;
shift_step=1;
projs=cryo_gen_projections(particle_size,n_particles,dummysnr,max_shift,shift_step);

% Put images in micrograph
for k=1:n_particles
    [i,j]=ind2sub([M M],idx(k));
    I=(i-1)*particle_size+1:i*particle_size;
    J=(j-1)*particle_size+1:j*particle_size;
    mgraph(I,J)=projs(:,:,k);
end
