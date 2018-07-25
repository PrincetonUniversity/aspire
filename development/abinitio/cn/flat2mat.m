function mat_arr = flat2mat(flat_arr,n)

[nr,nc,nd] = size(flat_arr);
assert((nr == 3) && (nc == 3) && nd == (nchoosek(n,2)));

mat_arr = zeros(3*n,3*n);

ind = 0;
for i=1:n
    for j=i+1:n
        ind = ind + 1;
        mat_arr((i-1)*3+1:i*3,(j-1)*3+1:j*3) = flat_arr(:,:,ind);
    end
end

end