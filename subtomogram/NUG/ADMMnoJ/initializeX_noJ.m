%{
INPUT:
    n: number of images
    d0,d1: vector containing dimensions of the representations
OUTPUT:
    X0: k-by-k block
    X1: (k+1)-by-(k+1) block
    Xq: quaternion block

Afonso Bandeira, Yutong Chen, Amit Singer
Princeton University, August 2016
%}
function [ X , Xq ] = initializeX_noJ( n , d )

%
% X0 = zeros(sum(d0.^2),round(n*(n+1)/2));
% 
% count = 0;
% for k = 1:length(d0)
%     for col = 1:d0(k)
%     for row = 1:d0(k)
%         
%         count = count + 1;
%         
%         if row == col
%             
%             X0(count,1:n) = ones(1,n);
%             
%             if k == 1
%                 X0(count,(n+1):end) = ones(1,round(n*(n-1)/2));
%             end
%             
%         end
%         
%     end
%     end
% end

%
% X1 = zeros(sum(d1.^2),round(n*(n+1)/2));
% 
% count = 0;
% for k = 1:length(d1)
%     for col = 1:d1(k)
%     for row = 1:d1(k)
%         
%         count = count + 1;
%         
%         if row == col
%             
%             X1(count,1:n) = ones(1,n);
%             
%             if k == 1
%                 X1(count,(n+1):end) = ones(1,round(n*(n-1)/2));
%             end
%             
%         end
%         
%     end
%     end
% end

% add
X = zeros(sum(d.^2),round(n*(n+1)/2));

count = 0;
for k = 1:length(d)
    for col = 1:d(k)
    for row = 1:d(k)
        
        count = count + 1;
        
        if row == col
            
            X(count,1:n) = ones(1,n);
            
            if k == 1
                X(count,(n+1):end) = ones(1,round(n*(n-1)/2));
            end
            
        end
        
    end
    end
end

%
Xq = zeros(16,round(n*(n-1)/2)); % 3*3 Y with additional 1 and zeros to 4*4
Xq(1,:) = ones(1,round(n*(n-1)/2)); % colunmized so each starts with 1

end









