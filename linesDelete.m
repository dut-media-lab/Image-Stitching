function [ inlier_line1, inlier_line2 ] = linesDelete( line1, line2, pts1, pts2 )
% outlier removal for line match
outlier_threshold = 3;%  if projection distance >3, outliers

% delete duplicate feature match
[h1, indh1] = sortrows(pts1');
[~,  indk1] = unique(h1,'rows');
matches1 = pts1(:, indh1(indk1));
matches2 = pts2(:, indh1(indk1));
[h2, indh2] = sortrows(matches2');
[~, indk2] = unique(h2, 'rows');
matches1 = matches1(:, indh2(indk2));
matches2 = matches2(:, indh2(indk2));

% outlier removal based on homography warp
init_H = calcHomo(matches1, matches2);
num_line = size(line1,1);
abc_line1 = [line1(:,4)-line1(:,2), line1(:,1)-line1(:,3), line1(:,3).*line1(:,2)-line1(:,1).*line1(:,4)];
abc_line2 = [line2(:,4)-line2(:,2), line2(:,1)-line2(:,3), line2(:,3).*line2(:,2)-line2(:,1).*line2(:,4)];

aux_line1 = [line1(:,1:2)';ones(1,num_line);line1(:,3:4)';ones(1,num_line)];
aux_line2 = [line2(:,1:2)';ones(1,num_line);line2(:,3:4)';ones(1,num_line)];

warp_p1 = init_H*aux_line1(1:3,:); warp_p1 = warp_p1./repmat(warp_p1(3,:),3,1);
warp_p2 = init_H*aux_line1(4:6,:); warp_p2 = warp_p2./repmat(warp_p2(3,:),3,1);

dist_p1 = abs(sum(warp_p1'.*abc_line2, 2))./sqrt(abc_line2(:,1).^2+abc_line2(:,2).^2);
dist_p2 = abs(sum(warp_p2'.*abc_line2, 2))./sqrt(abc_line2(:,1).^2+abc_line2(:,2).^2);
mean_dist = (dist_p1+dist_p2)./2;
inliers = mean_dist <= outlier_threshold;  

warp_p1_ = init_H\aux_line2(1:3,:); warp_p1_ = warp_p1_./repmat(warp_p1_(3,:),3,1);
warp_p2_ = init_H\aux_line2(4:6,:); warp_p2_ = warp_p2_./repmat(warp_p2_(3,:),3,1);

dist_p1_ = abs(sum(warp_p1_'.*abc_line1, 2))./sqrt(abc_line1(:,1).^2+abc_line1(:,2).^2);
dist_p2_ = abs(sum(warp_p2_'.*abc_line1, 2))./sqrt(abc_line1(:,1).^2+abc_line1(:,2).^2);
mean_dist_ = (dist_p1_+dist_p2_)./2;
inliers_ = mean_dist_ <= outlier_threshold;

line_inliers = inliers & inliers_;

inlier_line1 = line1(line_inliers,:);
inlier_line2 = line2(line_inliers,:);


end

