function [ lines_vs, lines_us, lines_ue ] = generateUV( img1, img2, init_H, theta, C1 ,C2 )
% given an orientation of img, rotate the original control vertices Mv
% obtain the orthogonal equidistant lines
% generating v-slope lines, u-slope lines, u-equal lines
% based on u sample points in non-overlapping region
multi = 2;
[M, N, ~] = size(img1);  % size of img
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];  % rotation matrix
off_center = [(N+1)/2; (M+1)/2] - R*[(N+1)/2; (M+1)/2];  % offset of rotated center r.t. original center旋转后图像中心的偏移
%% roughly calculate the overlapping boundary to filter u-sample points
warp_pts1 = init_H*[1;1;1]; warp_pts1 = warp_pts1(1:2)./warp_pts1(3);
left_or_right = warp_pts1(1)<=1;  [sz1, sz2, ~] = size(img2);
if left_or_right  % left is target, line x=1 is boundary
    inv_pts1 = init_H\[1;-sz1;1];  inv_pts1 = inv_pts1(1:2)./inv_pts1(3);
    inv_pts2 = init_H\[1;2*sz1;1]; inv_pts2 = inv_pts2(1:2)./inv_pts2(3);
else  % right is target, line x=N is boundary
    inv_pts1 = init_H\[sz2;-sz1;1];  inv_pts1 = inv_pts1(1:2)./inv_pts1(3);
    inv_pts2 = init_H\[sz2;2*sz1;1]; inv_pts2 = inv_pts2(1:2)./inv_pts2(3);
end
x1 = inv_pts1(1);  y1 = inv_pts1(2);  x2 = inv_pts2(1);  y2 = inv_pts2(2);

%% generating mesh grid (C1*C2) to optimize warped control vertices
[X, Y] = meshgrid(linspace(2-N, 2*N-1, 3*multi*(C2+1)-2), linspace(2-M, 2*M-1, 3*multi*(C1+1)-2)); % mesh grid index
lines_vs = zeros(2*size(X,2), size(X,1));
k=1;
% rotated vertical line  slope-preserving
for j=1:size(X,2)
    tmp_Mv = [X(:,j)'; Y(:,j)'];  % each vertical line
    tmp_line = R*tmp_Mv + repmat(off_center, 1,length(tmp_Mv)); % rotated vertical line
    inner_ind = (tmp_line(1,:)>=1) & (tmp_line(1,:)<=N) & (tmp_line(2,:)>=1) & (tmp_line(2,:)<=M); % rotated vertical line lies in img
    lines_vs(2*k-1:2*k, 1:sum(inner_ind)) = tmp_line(:,inner_ind);  % useful rotated vertical line
    lines_vs(2*k-1:2*k, end) = sum(inner_ind);  % number of sample points on the line
    k=k+1;
end 
% rotated horizontal line slope-preserving
lines_us = zeros(2*size(X,1), size(X,2));
k=1;
for i=1:size(X,1)
    tmp_Mv = [X(i,:); Y(i,:)];  % each horizontal line
    tmp_line = R*tmp_Mv + repmat(off_center, 1,length(tmp_Mv)); % rotated horizontal line
    inner_ind = (tmp_line(1,:)>=1) & (tmp_line(1,:)<=N) & (tmp_line(2,:)>=1) & (tmp_line(2,:)<=M); % rotated horizontal line lies in img
    lines_us(2*k-1:2*k, 1:sum(inner_ind)) = tmp_line(:,inner_ind);  % useful rotated horizontal line
    lines_us(2*k-1:2*k, end) = sum(inner_ind);  % number of sample points on the line
    k=k+1;
end 

lines_vs(all(sum(lines_vs,2),2)==0,:)=[];  % delete all-zero rows
lines_vs(:,all(sum(lines_vs,1),1)==0)=[];  % delete all-zero columns
lines_us(all(sum(lines_us,2),2)==0,:)=[];  % delete all-zero rows
lines_us(:,all(sum(lines_us,1),1)==0)=[];  % delete all-zero columns

%% filter u sample points (omit points in the overlapping region)
newlines_u = zeros(size(lines_us));
for i=1:2:size(lines_us,1)-1
    num_u = lines_us(i,end);
    vec_prod = (x1-lines_us(i,1:num_u)).*(y2-lines_us(i+1,1:num_u))+(lines_us(i+1,1:num_u)-y1).*(x2-lines_us(i,1:num_u));
    left = vec_prod>0;  right = vec_prod<0; %left = vec_prod>0;  right = vec_prod<0;
    index_prod = (left_or_right & left) | ((~left_or_right) & right);  % find the sample points in the non-overlapping region
    newlines_u(i:i+1, 1:sum(index_prod)) = lines_us(i:i+1, index_prod);
    newlines_u(i:i+1, end) = sum(index_prod); % number of sample points
end
lines_ue = newlines_u;  % u sample points to preserve equidistant
lines_ue(all(sum(lines_ue,2),2)==0,:)=[];  % delete all-zero rows
lines_ue(:,all(sum(lines_ue,1),1)==0)=[];  % delete all-zero columns

if isempty(lines_ue)
    for i=1:2:size(lines_us,1)-1
        num_u = lines_us(i,end);
        vec_prod = (x1-lines_us(i,1:num_u)).*(y2-lines_us(i+1,1:num_u))+(lines_us(i+1,1:num_u)-y1).*(x2-lines_us(i,1:num_u));
        left = vec_prod<0;  right = vec_prod>0; %left = vec_prod>0;  right = vec_prod<0;
        index_prod = (left_or_right & left) | ((~left_or_right) & right);  % find the sample points in the non-overlapping region
        newlines_u(i:i+1, 1:sum(index_prod)) = lines_us(i:i+1, index_prod);
        newlines_u(i:i+1, end) = sum(index_prod); % number of sample points
    end
    lines_ue = newlines_u;  % u sample points to preserve equidistant
    lines_ue(all(sum(lines_ue,2),2)==0,:)=[];  % delete all-zero rows
    lines_ue(:,all(sum(lines_ue,1),1)==0)=[];  % delete all-zero columns
end      
             
% figure,imshow(img1);
% hold on
% plot(lines_us(1:2:end-1,1:end-1), lines_us(2:2:end,1:end-1),'rx', 'LineWidth',2);
% plot(lines_ue(1:2:end-1,1:end-1), lines_ue(2:2:end,1:end-1),'gx', 'LineWidth',2);
% hold off
    
end

