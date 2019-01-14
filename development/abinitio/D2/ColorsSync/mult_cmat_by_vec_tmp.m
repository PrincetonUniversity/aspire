
%multiply color matrix by vector v on the fly
% function [out]=mult_cmat_by_vec(cperms,v,N)
%    cmat=@(x) apply_mat(cperms,x,N);
%    out=eigs(cmat,
% end

function [out]=mult_cmat_by_vec_pd(cperms,v,N)
    
    tperms=[1,2,3; 1,3,2; 2,1,3; 2,3,1; 3,1,2; 3,2,1]-1; %[1,2,3] permutations
    iperms=[1,2,3; 1,3,2; 2,1,3; 3,1,2; 2,3,1; 3,2,1]-1; %inverse permutations
    out=zeros(length(v),1);
    p_n=zeros(3,1);
    for i=1:N-2
        for j=i+1:N-1
            for k=j+1:N
                ijk=trip_idx(i,j,k,N);
                ij=3*uppertri_ijtoind(i,j,N)-2;
                ik=3*uppertri_ijtoind(i,k,N)-2;
                jk=3*uppertri_ijtoind(j,k,N)-2;
                %extract permutation indexes from cperms
                n=cperms(ijk);
                p_n(1)=floor(n/100);
                p_n(3)=mod(n,10);
                p_n(2)=(n-p_n(1)*100-p_n(3))/10;
                
                %multiply vector by color matrix
                %Upper triangular part
                p=tperms(p_n(1),:)+ik;
                out(ij)=out(ij)-v(p(2))-v(p(3))+v(p(1));
                out(ij+1)=out(ij+1)-v(p(1))-v(p(3))+v(p(2));
                out(ij+2)=out(ij+2)-v(p(1))-v(p(2))+v(p(3));
                p=tperms(p_n(2),:)+jk;
                out(ij)=out(ij)-v(p(2))-v(p(3))+v(p(1));
                out(ij+1)=out(ij+1)-v(p(1))-v(p(3))+v(p(2));
                out(ij+2)=out(ij+2)-v(p(1))-v(p(2))+v(p(3));
                p=iperms(p_n(3),:)+jk;
                out(ik)=out(ik)-v(p(2))-v(p(3))+v(p(1));
                out(ik+1)=out(ik+1)-v(p(1))-v(p(3))+v(p(2));
                out(ik+2)=out(ik+2)-v(p(1))-v(p(2))+v(p(3));
                %Lower triangular part
                p=iperms(p_n(1),:)+ij;
                out(ik)=out(ik)-v(p(2))-v(p(3))+v(p(1));
                out(ik+1)=out(ik+1)-v(p(1))-v(p(3))+v(p(2));
                out(ik+2)=out(ik+2)-v(p(1))-v(p(2))+v(p(3));
                p=iperms(p_n(2),:)+ij;
                out(jk)=out(jk)-v(p(2))-v(p(3))+v(p(1));
                out(jk+1)=out(jk+1)-v(p(1))-v(p(3))+v(p(2));
                out(jk+2)=out(jk+2)-v(p(1))-v(p(2))+v(p(3));
                p=tperms(p_n(3),:)+ik;
                out(jk)=out(jk)-v(p(2))-v(p(3))+v(p(1));
                out(jk+1)=out(jk+1)-v(p(1))-v(p(3))+v(p(2));
                out(jk+2)=out(jk+2)-v(p(1))-v(p(2))+v(p(3));
                
            end
        end
    end
end




