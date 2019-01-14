
xy_proj=squeeze(bad_pairs(1:2,:,:));
norms=sqrt(sum(xy_proj.^2,1));
xy_proj=xy_proj./repmat(norms,2,1);
tmp=squeeze(dot(xy_proj(:,1,:),xy_proj(:,2,:),1));

xz_proj=squeeze(bad_pairs([1,3],:,:));
norms=sqrt(sum(xz_proj.^2,1));
xz_proj=xz_proj./repmat(norms,2,1);
tmp=[tmp,squeeze(dot(xz_proj(:,1,:),xz_proj(:,2,:),1))];

yz_proj=squeeze(bad_pairs(1:2,:,:));
norms=sqrt(sum(yz_proj.^2,1));
yz_proj=yz_proj./repmat(norms,2,1);
tmp=[tmp,squeeze(dot(yz_proj(:,1,:),yz_proj(:,2,:),1))];


tmp_err=squeeze(sum(reshape((tmp-tmp2).^2,9,size(tmp,3)),1));
tmp235=reshape((tmp-tmp2).^2,9,size(tmp,3));
