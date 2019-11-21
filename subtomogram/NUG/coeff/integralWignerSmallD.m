function H = integralWignerSmallD(t)

H = cell(t+1,1);

for kk = 1:(t+1)
    
    k = kk - 1;
    
    wk = wignerd(k,pi/2);
    
    mat = zeros(2*k+1);
    
    for m1 = (-k):k
    for m2 = (-k):k
        
        row = m1 + (k+1); col = m2 + (k+1);
        
        val = 0;
        
        for l = (-k):k
            
            if l == 0
                val = val + 2*wk(k+1,row)*wk(k+1,col);
            elseif l == 1
                val = val - 1i*pi/2*wk(k+2,row)*wk(k+2,col);
            elseif l == -1
                val = val + 1i*pi/2*wk(k,row)*wk(k,col);
            else
                val = val + ((1+exp(-1i*l*pi))/(1-l^2))*wk(l+(k+1),row)*wk(l+(k+1),col);
            end
            
        end
        
        mat(row,col) = (-1)^(m1+m2)*1i^(m1-m2)*val;
        
    end
    end
    
    H{kk} = real(mat);
    
end

end









