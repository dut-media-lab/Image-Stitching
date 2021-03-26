function [ h, A, T1, T2 ] = calcHomoPointLine( pts1, pts2, line1, line2 )
% estimate homography based on feature matches and line segments matches
% pts1, pts2: 2*N, N matches (x1,x2,...;y1,y2,...)
% line1, line2: M*4,  M matches (x1,y1,x2,y2;...) two endpoints of line
% segments
num_pts = size(pts1,2);  num_line = size(line1,1);

%% point-centric normalization before SVD r.t Dual-feature
% for target image
re_line1 = reshape(line1',2,2*num_line);
[ norm_pts_line1, T1 ] = normalise2dpts([[pts1,re_line1]; ones(1,num_pts+2*num_line)]);
norm_pts1 = norm_pts_line1(1:2, 1:num_pts);
norm_line1 = reshape(norm_pts_line1(1:2, num_pts+1:end),4,num_line);

% for reference image
abc_line2 = [line2(:,4)-line2(:,2), line2(:,1)-line2(:,3),...
                line2(:,3).*line2(:,2)-line2(:,1).*line2(:,4)]; % line's function: ax+by+c=0
tmp_line2 = abc_line2(:,1:3)./repmat(abc_line2(:,3),1,3);
%[ norm_lines2, T3 ]=normaliselines(tmp_line2);
% least-square method
x0 = [1,0,0]; options.Algorithm = 'levenberg-marquardt';
x = lsqnonlin(@(x)myfun(x, pts2, tmp_line2), x0,[],[],options);%
T2 = [x(1),0,x(2);0,x(1),x(3);0,0,1];
norm_pts2 = T2*[pts2;ones(1,num_pts)];
norm_line2 = [x(1).*abc_line2(:,1),x(1).*abc_line2(:,2), x(1)^2.*abc_line2(:,3)-x(1)*x(2).*abc_line2(:,1)-x(1)*x(3).*abc_line2(:,2)];
            
%% generating feature matches value matrix A (||Ah||=0)  
normA_pts = zeros(2*num_pts,9);
normB_line = zeros(2*num_line,9);

x1 = norm_pts1(1,:)'; y1 = norm_pts1(2,:)';
x2 = norm_pts2(1,:)'; y2 = norm_pts2(2,:)';
normA_pts(1:2:2*num_pts-1,1:3) = [x1, y1, ones(num_pts,1)];
normA_pts(1:2:2*num_pts-1,7:9) = [-x2.*x1, -x2.*y1, -x2.*ones(num_pts,1)];
normA_pts(2:2:2*num_pts,4:6)   = [x1, y1, ones(num_pts,1)];
normA_pts(2:2:2*num_pts,7:9)   = [-y2.*x1, -y2.*y1, -y2.*ones(num_pts,1)];

%% generating line segments matches value matrix B (||Bh||=0) 
u0 = norm_line1(1,:)';  v0=norm_line1(2,:)';
u1 = norm_line1(3,:)';  v1=norm_line1(4,:)';
a2=norm_line2(:,1); b2=norm_line2(:,2); c2=norm_line2(:,3);
k=1./(a2.^2+b2.^2);
normB_line(1:2:2*num_line-1,:)=repmat(k,1,9).*[a2.*u0,a2.*v0,a2,b2.*u0,b2.*v0,b2,c2.*u0,c2.*v0,c2];
normB_line(2:2:2*num_line,:)  =repmat(k,1,9).*[a2.*u1,a2.*v1,a2,b2.*u1,b2.*v1,b2,c2.*u1,c2.*v1,c2];

%%  Solve homography with ||[A;B]h=0|| 
norm_C = [normA_pts; normB_line];  
[~,~,V] = svd(norm_C,0);
h = reshape(V(:,9),3,3)';

A = norm_C;


end
 
function F = myfun(x, pts2, tmp_line2) % least-square cost function

a=tmp_line2(:,1); b=tmp_line2(:,2); c=tmp_line2(:,3);

F(1) = mean(sqrt((x(1)*pts2(1,:)+x(2)).^2+(x(1)*pts2(2,:)+x(3)).^2))-sqrt(2); % "average point [1,1,1]"
F(2) = mean(abs(x(1)*x(1).*c-x(1)*x(2).*a-x(1)*x(3).*b)./sqrt((x(1).*a).^2+(x(1).*b).^2))-sqrt(1/2);

end

