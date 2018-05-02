function [projs,refq] = remove_eq_images(projs,refq)

log_message('removing equator images');
nImages = size(refq,2);
assert(nImages == size(projs,3));

is_eq = zeros(1,nImages);
for i=1:nImages
    
   rot = q_to_rot(refq(:,i)).';
   if ( (abs (acosd(rot(3,3))) - 90) < 10 )
       is_eq(i) = 1;
   end
end

eq_inds = find(is_eq == 1);

refq(:,eq_inds) = [];
projs(:,:,eq_inds) = [];

end