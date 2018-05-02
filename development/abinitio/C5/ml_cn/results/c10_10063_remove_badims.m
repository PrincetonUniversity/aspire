vol = ReadMRC('/home/yoel/scratch/10063/particles/ring10.mrcs');

bad_inds = [];
for i=1:size(vol,3)
    voli = vol(:,:,i);
    is_all_nonzero_col = all(any(voli,1));
    is_all_nonzero_row = all(any(voli,2));

    if ~is_all_nonzero_col || ~is_all_nonzero_row || norm(voli(:)) < 0.1
        bad_inds(end+1) = i;
    end

end

vol(:,:,bad_inds) = [];
WriteMRC(vol,1,'datasets/10063_c10/ring10_wo_blackims.mrcs');
