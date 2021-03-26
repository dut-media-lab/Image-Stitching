function [dist, H] = homography_res(H, X)
    
    H = reshape(H,3,3);
   
    
    x1 = X(1:3,:);   % Extract x1 and x2 from x
    x2 = X(4:6,:);    
    n = size(x1,2);
    
    % Calculate, in both directions, the transfered points    
    Hx1    = H*x1;
    invHx2 = H\x2;
    
    % Normalise so that the homogeneous scale parameter for all coordinates
    % is 1.
    
    x1     = hnormalise(x1);
    x2     = hnormalise(x2);     
    Hx1    = hnormalise(Hx1);
    invHx2 = hnormalise(invHx2); 
    
    dist = sum((x1-invHx2).^2)  + sum((x2-Hx1).^2);
    
    dist = reshape(dist,n,1);
    
    H = H(:);

end



% function [dist, P] = homography_res(P, X)
% % Function to calculate distances between a homography and a set of point
% % correspondences. Implmenents the symmetric transfer error. See
% % [Harley and Zisserman 2nd ed., pg 94].
% % X is a 4 by npts matrix.
% % X(1,i) is the x coordinate in the 1st image of the i-th correspondence.
% % X(2,i) is the y coordinate in the 1st image of the i-th correspondence.
% % X(3,i) is the x coordinate in the 2nd image of the i-th correspondence.
% % X(4,i) is the y coordinate in the 2nd image of the i-th correspondence.
% 
% [dummy npts] = size(X);
% dist = zeros(npts,1);    
% 
% % Forward.
% H = reshape(P,3,3)';
% hom = H*[ X(3:4,:) ; ones(1,size(X,2)) ];
% inhom = [ hom(1,:)./hom(3,:) ; hom(2,:)./hom(3,:) ];
% for i=1:npts
% dist(i) = norm(X(1:2,i) - inhom(:,i))^2;
% end
% 
% % Reverse.
% % invH = inv(H);
% % hom = invH*[ X(1:2,:) ; ones(1,size(X,2)) ];
% 
% hom = H\[ X(1:2,:) ; ones(1,size(X,2)) ];
% 
% inhom = [ hom(1,:)./hom(3,:) ; hom(2,:)./hom(3,:) ];
% for i=1:npts
%     dist(i) = dist(i) + norm(X(3:4,i) - inhom(:,i))^2;
% end
%     
% end
