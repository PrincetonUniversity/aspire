function projs = get_ims_from_stack(mrc_stack,im_indeces)

%TODO: get the stack size (i.e. the size of each each image along with the
%number of images) directly from the mrc
% TODO: assert that no input im index is bigger than the number of images
% in th mrc file
stack = imagestackReader(mrc_stack);
tmp_img = stack.getImage(1);
sz = size(tmp_img,1);

projs = zeros(sz,sz,numel(im_indeces));
msg = [];
for k=1:nImages
    t1 = clock;
    ind = im_indeces(k);
    projs(:,:,k) = stack.getImage(ind);
       
    %%%%%%%%%%%%%%%%%%% debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t2 = clock;
    t = etime(t2,t1);
    bs = char(repmat(8,1,numel(msg)));
    fprintf('%s',bs);
    msg = sprintf('k=%3d/%3d  t=%7.5f',k,nImages,t);
    fprintf('%s',msg);
    %%%%%%%%%%%%%%%%%%% end of debug code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

end