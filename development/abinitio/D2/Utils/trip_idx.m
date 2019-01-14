
function res=trip_idx(i,j,k,N)    
    %res=0.25*nchoosek(N-i,3)-0.5*(j-2)*(2*(N-i)+j-3)-(k-j);
%     res=nchoosek(N,3)-nchoosek(N-i,3)-nchoosek(N-j+1,2)+k-j;
    res=N*(N-1)*(N-2)/6-(N-i-2).*(N-i-1).*(N-i)/6-(N-j).*(N-j+1)/2+k-j;
 %   res=(N-1)*(N-2)*(N-3)/6-(N-i-4)*(N-i-3)*(N-i-2)/6-(N-j-2)*(N-j-1)/2+k-j-1;Wrong
   % res=N*(N-1)*(N-2)/6-(N-i-3)*(N-i-2)*(N-i-1)/6-(N-j-1)*(N-j)/2+k-j-1;
end