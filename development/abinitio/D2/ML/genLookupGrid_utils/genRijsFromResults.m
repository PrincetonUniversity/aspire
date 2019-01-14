
%generate Rij pairs from cryo_clmatrix results
function [Rijs]=genRijsFromResults(cl_data,cl_idx)

n=cl_data.n;
n2=cl_data.n2;
ntheta=cl_data.ntheta;
dtheta=2*pi/ntheta;
RiRj0_lookup=cl_data.RiRj0;
npairs=nchoosek(n,2);
npairs2=n*n2;
nproj_pairs=length(cl_idx);

[thetas_j,thetas_i,p_idx]=...
    ind2sub([ntheta,ntheta,2*npairs+2*npairs2],cl_idx);
% [~,~,ij_partition]=histcounts(p_idx,...
%     [1,npairs+1,2*npairs+1,2*npairs+npairs2+1,2*npairs+2*npairs2]);

%DEBUG CODE
p_idx_deb=ceil(cl_idx/ntheta^2); %index of pair from all pairs
[thetas_j_deb,thetas_i_deb]=...
    ind2sub([ntheta,ntheta],cl_idx-(p_idx_deb-1)*ntheta^2);
%END
    
%Generate Rijs from results
Ris0=squeeze(RiRj0_lookup(:,:,1,p_idx));
Rjs0=squeeze(RiRj0_lookup(:,:,2,p_idx));
cosines_i=cos(2*pi-(dtheta*(thetas_i-1)));
cosines_j=cos(2*pi-(dtheta*(thetas_j-1)));
sines_i=sin(2*pi-(dtheta*(thetas_i-1)));
sines_j=sin(2*pi-(dtheta*(thetas_j-1)));
inplane_rots_i=zeros(3,3,nproj_pairs);
inplane_rots_j=zeros(3,3,nproj_pairs);

inplane_rots_i(1,1,:)=cosines_i;
inplane_rots_i(2,2,:)=cosines_i;
inplane_rots_i(1,2,:)=-sines_i;
inplane_rots_i(2,1,:)=sines_i;
inplane_rots_i(3,3,:)=1;

inplane_rots_j(1,1,:)=cosines_j;
inplane_rots_j(2,2,:)=cosines_j;
inplane_rots_j(1,2,:)=-sines_j;
inplane_rots_j(2,1,:)=sines_j;
inplane_rots_j(3,3,:)=1;

Ris=multiprod(Ris0,inplane_rots_i);
Rjs=multiprod(Rjs0,inplane_rots_j);
Ris_t=multi_transpose(Ris);

clearvars Ris0 Rjs0 Ris inplane_rots_i inplane_rots_j cosines_i...
          cosines_j sines_i sines_j thetas_j thetas_i p_idx

Rijs=multiprod(Ris_t,Rjs);
gxRjs=Rjs;
gxRjs(2:3,:)=-Rjs(2:3,:);
Rijs_x=multiprod(Ris_t,gxRjs);
gyRjs=Rjs;
gyRjs([1,3],:)=-Rjs([1,3],:);
Rijs_y=multiprod(Ris_t,gyRjs);
gzRjs=Rjs;
gzRjs(1:2,:)=-Rjs(1:2,:);
Rijs_z=multiprod(Ris_t,gzRjs);

Rijs=cat(4,Rijs,Rijs_x,Rijs_y,Rijs_z);
Rijs=permute(Rijs,[1,2,4,3]);
return;
%% DEBUG CODE- Test Rijs vs common lines
cls=zeros(4,2,nproj_pairs);
for l=1:4
    cls(l,1,:)=mod(atan2(Rijs(1,3,l,:),-Rijs(2,3,l,:))+2*pi,2*pi);
    cls(l,2,:)=mod(atan2(-Rijs(3,1,l,:),Rijs(3,2,l,:))+2*pi,2*pi);
end

% for i=1:nproj_pairs
%     for l=1:4
%         cls(l,1,ind)=mod(atan2(Rijs(1,3,l,ind),-Rijs(2,3,l,ind))+2*pi,2*pi);
%         cls(l,2,ind)=mod(atan2(-Rijs(3,1,l,ind),Rijs(3,2,l,ind))+2*pi,2*pi);
%     end
% end

cls=reshape(cls,8,nproj_pairs);
cls=mod(round(cls*180/pi),360);
load('thetas_500_10deg.mat');
if exist('thetas','var')
    s=size(thetas);
    thetas=reshape(thetas,4,2,s(3)*s(4));
    thetas=reshape(thetas,8,size(thetas,3));
    diff=cls-thetas;
    diff=sum(abs(diff(:)));
    disp(['diff==',num2str(diff)]);
end





    