
%Get random rotations from lookup grid
%Input: ns = number of rotations to sample
%       nips = number of inplane rotations for each sample
function [test_data]=testRotFromGrid(lookup_data,ns,nips)
irg=lookup_data.inplane_rotated_grid;
dims=size(irg);
nr=dims(4);
nip=dims(3);

if ns>nr
    error('Number of samples is larger than the grid');
end
if nips>nip
    error('Number of in-plane samples is larger than number of inplane rotations in the grid');
end
sampledGrid=zeros(3,3,nips,ns);

s_idx=randperm(nr,ns);
s_idx=sort(s_idx);
ip_idx=zeros(ns,nips);
for k=1:ns
    ip_idx(k,:)=randperm(nip,nips);
    sampledGrid(:,:,:,k)=irg(:,:,ip_idx(k,:),s_idx(k));
end
test_data.rots_grid=squeeze(sampledGrid);
test_data.s_idx=s_idx';
test_data.ip_idx=ip_idx;



