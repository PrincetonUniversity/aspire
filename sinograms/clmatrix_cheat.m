function [clmatrix,clcorr,cltheta]=clmatrix_cheat(rots,n_theta)
%
% Build common lines matrix using the true rotations corresponding to the
% projections orientations. Each projection has n_theta rays.
% clcorr is set to 1.0e-8.
%
% clmatrix contains the indices of the common lines between each pair of
% projections. cltheta contains the exact angle of the common lines and
% thus does not contain any discretization errors.
%
% Yoel Shkolnisky, October 2008.
%
% Revised, Y.S. June 2014.


N=size(rots,3);

clmatrix=zeros(N);   % common lines matrix
clcorr=zeros(N);     % correlation coefficient for ach common line
cltheta=zeros(N);    % angles of common line pairs

for k1=1:N-1
    R1=rots(:,:,k1);
    R1=R1.';
    
    % Rotated coordinate system is the columns of the rotation matrix.
    % We use explicit multiplication to show which column corresponds to x,
    % which to y, and which to z. See commonline_euler for the difference.
    X1=R1(:,1);
    Y1=R1(:,2);
    Z1=R1(:,3);
    
    XY1=[X1 Y1];
    
    for k2=k1+1:N
        %%%%%%%%%%%%%%%%%%%%%%%
                
        R2=rots(:,:,k2);
        R2=R2.';
        
        X2=R2(:,1);
        Y2=R2(:,2);
        Z2=R2(:,3);
        
        Z3=[Z1(2)*Z2(3)-Z1(3)*Z2(2);...
            Z1(3)*Z2(1)-Z1(1)*Z2(3);...
            Z1(1)*Z2(2)-Z1(2)*Z2(1)];
        
        % Make sure the projections are not too close.
        if norm(Z3)<1.0e-8
            warning('GCAR:normTooSmall','Images have same orientation');
        end
        
        Z3=Z3./norm(Z3);
        
        % Compute coordinates of the common-line in each local coordinate system.
        XY2=[X2 Y2];
        c1=(Z3.')*XY1;
        c2=(Z3.')*XY2;
        
        % Verify that the common-line is indeed common to both planes. The
        % following warning should never happen! Just to make sure nothing went
        % terribly wrong.
        ev1=XY1*c1(:)-Z3;
        ev2=XY2*c2(:)-Z3;
        
        if (norm(ev1)/norm(Z3)>1.0e-12) || (norm(ev2)/norm(Z3)>1.0e-12)
            warning('GCAR:largeErrors',...
                'Common line is not common. Error1 = %e, Error2 = %e',...
                norm(ev1)/norm(Z3),norm(ev2)/norm(Z3));
        end
        
        % Compute angle of the common line at each projection's coordinate system
        theta1=atan2(c1(2),c1(1));
        theta2=atan2(c2(2),c2(1));
        
        PI=4*atan(1.0);
        theta1=theta1+PI; % Shift from [-pi,pi] to [0,2*pi].
        theta2=theta2+PI;
        
        idx1=theta1/(2*PI)*n_theta;
        idx2=theta2/(2*PI)*n_theta;
        
        idx1=mod(round(idx1),n_theta);
        idx2=mod(round(idx2),n_theta);
        %%%%%%%%%%%%%%%%%%%%%%%%
        clmatrix(k1,k2)=idx1+1;
        clmatrix(k2,k1)=idx2+1;
        clcorr(k1,k2)=1.0e-8;
        cltheta(k1,k2)=theta1;
        cltheta(k2,k1)=theta2;
    end
end
