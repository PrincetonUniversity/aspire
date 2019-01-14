function [Jified_rot]=multi_Jify(in)
    
    dims=size(in);
    l=length(dims);
    Jified_rot=reshape(in,[9 dims(3:l)]);
    Jified_rot([3 6 7 8],:,:)=-Jified_rot([3 6 7 8],:,:);
    Jified_rot=reshape(Jified_rot,[3 3 dims(3:l)]);

end