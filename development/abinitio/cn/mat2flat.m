function flat_arr = mat2flat(mat_arr,n)

[nr,nc] = size(mat_arr);
assert( (n*3 == nr) && (n*3 == nc) );

flat_arr = zeros(3,3,nchoosek(n,2));

ind = 0;
for i=1:n
    for j=i+1:n
        ind = ind + 1;
        flat_arr(:,:,ind) = mat_arr((i-1)*3+1:i*3,(j-1)*3+1:j*3);
    end
end

end