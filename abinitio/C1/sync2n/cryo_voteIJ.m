function [goodK3,peakh,alpha]=cryo_voteIJ(clmatrix,L,k1,k2,K3,rots_ref,is_perturbed)
%
% Apply the voting algorithm for images k1 and k2. clmatrix is the common
% lines matrix, constructed using angular reslution L. K3 are the images to
% use for voting of the pair (k1,k2).
% verbose  1 for detailed massages, 0 for silent execution (default).
%
% The function returns goodK3, the list of all third projections in the
% peak of the histogram corresponding to the pair (k1,k2), and the
% estimated angle alpha between them.
%
% Yoel Shkolnisky, September 2010.
%
% Revisions:
% Y.S. March 01, 2015  Optimize the code by replacing the loop on k3 with
%       vectorized operations. Also, optimize histogram smoothing by
%       replacing the loop with vectorized opeartions.


%K=size(clmatrix,1);

% Parameters used to compute the smoothed nagle histogram. Same as in 
% histfilter_v6 (for compatibility).
ntics=60;  
x=linspace(0,180,ntics);
%h=zeros(size(x));

% %% START reference code
% l_idx=zeros(3); % The matrix of common lines between the triplet (k1,k2,k3).
% l_idx(1,2) = clmatrix(k1,k2);
% l_idx(2,1) = clmatrix(k2,k1);
% 
% phis=zeros(numel(K3),2); % Angle between k1 and k2 induced by each third
%                    % projection k3. The first column is the cosine of
%                    % that angle, The second column is the index k3 of
%                    % the projection that creates that angle.
% 
% idx=1;
% 
% rejected=zeros(numel(K3),1);
% rejidx=0;
% 
% for k3=K3
%     if (k1~=k2) && (k1~=k3) && clmatrix(k1,k2)~=0 && clmatrix(k1,k3)~=0 && clmatrix(k2,k3)~=0
%         % some of the entries in clmatrix may be zero if we cleared
%         % them due to small correlation, or if for each image
%         % we compute intersections with only some of the other
%         % images.
%         % l_idx=clmatrix([k1 k2 k3],[k1 k2 k3]);
%         %
%         % Note that as long as the diagonal of the common lines matrix is
%         % zero, the conditions (k1~=k2) && (k1~=k3) are not needed, since
%         % if k1==k2 then clmatrix(k1,k2)==0 and similarly for k1==k3 or
%         % k2==k3. Thus, the previous voting code (from the JSB paper) is
%         % correct even though it seems that we should test also that
%         % (k1~=k2) && (k1~=k3) && (k2~=k3), and only (k1~=k2) && (k1~=k3)
%         % as tested there. 
%         
%         l_idx(1,3)=clmatrix(k1,k3);
%         l_idx(3,1)=clmatrix(k3,k1);
%         l_idx(2,3)=clmatrix(k2,k3);
%         l_idx(3,2)=clmatrix(k3,k2);
%         
%         % theta1 is the angle on C1 created by its intersection with C3 and C2.
%         % theta2 is the angle on C2 created by its intersection with C1 and C3.
%         % theta3 is the angle on C3 created by its intersection with C2 and C1.
%         theta1 = (l_idx(1,3) - l_idx(1,2))*2*pi/L;
%         theta2 = (l_idx(2,1) - l_idx(2,3))*2*pi/L;
%         theta3 = (l_idx(3,2) - l_idx(3,1))*2*pi/L;
%         
%         c1=cos(theta1);
%         c2=cos(theta2);
%         c3=cos(theta3);
%         % Each common-line corresponds to a point on the unit sphere. Denote the
%         % coordinates of these points by (Pix, Piy Piz), and put them in the matrix
%         %   M=[ P1x  P2x  P3x ; ...
%         %       P1y  P2y  P3y ; ...
%         %       P1z  P2z  P3z ].
%         %
%         % Then the matrix
%         %   C=[ 1 c1 c2 ;...
%         %       c1 1 c3 ;...
%         %       c2 c3 1],
%         % where c1,c2,c3 are given above, is given by C=M.'*M.
%         % For the points P1,P2, and P3 to form a triangle on the unit shpere, a
%         % necessary and sufficient condition is for C to be positive definite. This
%         % is equivalent to
%         %      1+2*c1*c2*c3-(c1^2+c2^2+c3^2)>0.
%         % However, this may result in a traingle that is too flat, that is, the
%         % angle between the projections is very close to zero. We therefore use the
%         % condition below
%         %       1+2*c1*c2*c3-(c1^2+c2^2+c3^2) > 1.0e-5
%         % This ensures that the smallest singular value (which is actually
%         % controlled by the determinant of C) is big enough, so the matrix is far
%         % from singular. This condition is equivalent to computing the singular
%         % values of C, followed by checking that the smallest one is big enough.
%         
%         if 1+2*c1*c2*c3-(c1^2+c2^2+c3^2) > 1.0e-5
%             
%             cos_phi2 = (c3-c1*c2)/(sin(theta1)*sin(theta2));
%             if abs(cos_phi2)>1
%                 if abs(cos_phi2)-1>1.0e-12
%                     warning('GCAR:numericalProblem','cos_phi2>1. diff=%7.5e...Setting to 1.',abs(cos_phi2)-1);
%                 else
%                     cos_phi2=sign(cos_phi2);
%                 end
%             end
%             
%             phis(idx,1)=cos_phi2;
%             phis(idx,2)=k3;
%             idx=idx+1;
%         else
%             rejidx=rejidx+1;
%             rejected(rejidx)=k3;
%         end
%     end
% end
% 
% phis=phis(1:idx-1,:);
% rejected=rejected(1:rejidx);
% %% END reference code

% Optimized implementation of the above reference code. See comments above
% for details.
phis=zeros(numel(K3),2);
rejected=zeros(numel(K3),1);
idx=1;
rejidx=1;

if (k1~=k2) && clmatrix(k1,k2)~=0
    l_idx12 = clmatrix(k1,k2);
    l_idx21 = clmatrix(k2,k1);

    K3=K3((k1~=K3) & (clmatrix(k1,K3)~=0) & (clmatrix(k2,K3)~=0));
    
    l_idx13=clmatrix(k1,K3);
    l_idx31=clmatrix(K3,k1);
    l_idx23=clmatrix(k2,K3);
    l_idx32=clmatrix(K3,k2);
    
    theta1 = (l_idx13 - l_idx12)*2*pi/L;
    theta2 = (l_idx21 - l_idx23)*2*pi/L;
    theta3 = (l_idx32 - l_idx31)*2*pi/L;
    
    theta1=theta1(:);
    theta2=theta2(:);
    theta3=theta3(:);
    
    c1=cos(theta1);
    c2=cos(theta2);
    c3=cos(theta3);
    
    cond=1+2.*c1.*c2.*c3-(c1.^2+c2.^2+c3.^2);
    
    goodidx=find(cond>1.0e-5);
    badidx=find(cond<=1.0e-5);
    
    cos_phi2 = (c3(goodidx)-c1(goodidx).*c2(goodidx))./(sin(theta1(goodidx)).*sin(theta2(goodidx)));
    checkidx=find(abs(cos_phi2)>1);
    if any(abs(cos_phi2)-1>1.0e-12)
        warning('GCAR:numericalProblem','cos_phi2>1. diff=%7.5e...Setting to 1.',abs(cos_phi2)-1);
    elseif ~isempty(checkidx)
        cos_phi2(checkidx)=sign(cos_phi2(checkidx));
    end
    
    
    phis(idx:idx+numel(goodidx)-1,1)=cos_phi2;
    phis(idx:idx+numel(goodidx)-1,2)=K3(goodidx);
    idx=idx+numel(goodidx);
    
    rejected(rejidx:rejidx+numel(badidx)-1)=K3(badidx);
    rejidx=rejidx+numel(badidx);
    
end

phis=phis(1:idx-1,:);
rejected=rejected(1:rejidx-1);

% End of optimized implementation

goodK3=[];
peakh=-1;
alpha=-1;

if idx>1
    % Compute the histogram of the angles between projections k1
    % and k2.
    
    angles=acos(phis(:,1))*180/pi;
    % Angles that are up to 10 degrees apart are considered
    % similar. This sigma ensures that the width of the density
    % estimation kernel is roughly 10 degrees. For 15 degress, the
    % value of the kernel is negligible.
    
    %sigma=2.64;
    sigma=3.0; % For compatibility with histfilter_v6

%%% START reference code
%   tic;
%     for j=1:numel(x)
%         h(j)=sum(exp(-(x(j)-angles).^2/(2*sigma.^2)));
%     end
%    Told=toc;
%%% END reference code

% START optimization of the reference code above
%     tic;
    h=2*angles*x;
    h=bsxfun(@minus,h,angles.^2);
    h=bsxfun(@minus,h,x.^2);
    h=exp(h./(2*sigma.^2));    
    h=sum(h);
%     Tnew=toc;
%     if norm(h2-h)/norm(h)>1.0e-12
%         aaa=1;
%     end
%     fprintf('Told/Tnew=%5.3f\n',Told/Tnew);
    
% END optimization

    % We assume that at the location of the peak we get the true angle
    % between images k1 and k2. Find all third images k3, that induce an
    % angle between k1 and k2 that is at most 10 off the true
    % angle. Even for debugging, don't put a value that is smaller than two
    % tics, since the peak might move a little bit due to wrong k3 images
    % that accidentally fall near the peak.
    [peakh,peakidx]=max(h);
    idx=find(abs(angles-x(peakidx))<360/ntics);
    goodK3=phis(idx,2);
    alpha=phis(idx,1);
    
%     plot(x,h)
%     hold on;    
%     clr=repmat([1 0 0],numel(angles),1);   
%     clr(idx,:)=repmat([0 1 0],numel(idx),1);
%     scatter(angles,zeros(size(angles)),20,clr);
%     hold off;


    if ~isscalar(rots_ref)
        % For debugging, compute the true angle between the images, and compare to
        % the estimated one. The following debugging code is correct only when
        % p=1, for otherwise, it is hard to predict how the errors would affect
        % the results.
        R1=rots_ref(:,:,k1)';
        R2=rots_ref(:,:,k2)';
        alpharef=dot(R1(:,3),R2(:,3));


        if ~isscalar(is_perturbed)
            if ~is_perturbed(k1,k2) % Check the angle only for correct common lines.
                % Otherwise, we expect grabageanyway.
                if (max(abs(acosd(alpha-alpharef)-90))>360/ntics) && ...
                        (max(abs(acosd(alpha+alpharef)-90))>360/ntics)
                    warning('Voted the wrong angle');
                end

                % Check that we have found the correct third images, and only those.
                % Note that it is possible for wrong k3 images to appear in the list
                % due to an angle that is mistakenly correct. No warnings
                % should appear if there are no errors in the common lines
                % matrix.
                % If the pair (k1,k2) is correct, we are looking for all images k3
                % such that (k1,k3) and (k2,k3) are correct. from this list
                % we subtract the list of images rejected due to too small
                % triangle. The resulting list should be identical to the list
                % goodK3.

                goodcount=0;
                goodK3ref=zeros(numel(K3),1);

                for k3=K3
                    if (k1~=k2) && (k2~=k3) && (k1~=k3)
                        kk1=min(k1,k3);
                        kk2=max(k1,k3);
                        kk3=min(k2,k3);
                        kk4=max(k2,k3);

                        if kk2>kk1 && kk4>kk3
                            if ~is_perturbed(kk1,kk2) && ~is_perturbed(kk3,kk4)
                                goodcount=goodcount+1;
                                goodK3ref(goodcount)=k3;
                            end
                        end
                    end
                end
                goodK3ref=goodK3ref(1:goodcount);
                goodK3ref=setdiff(goodK3ref,rejected);
                diff=setxor(goodK3ref,goodK3);
                if ~isempty(diff)
                    warning('Voted the wrong images');
                end
            end
        end
    end
end

        