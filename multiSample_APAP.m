function [matches_1, matches_2] = multiSample_APAP(pts1, pts2, img1, img2)   

global fitfn resfn degenfn psize numpar
fitfn = 'homography_fit';
resfn = 'homography_res';
degenfn = 'homography_degen';
psize = 4;
numpar = 9;

% normalise point distribution.
data_orig = [ pts1(1:2,:) ; ones(1,size(pts1,2)) ; pts2(1:2,:) ; ones(1,size(pts2,2)) ];
[data_norm_img1, ~] = normalise2dpts(data_orig(1:3,:));
[data_norm_img2, ~] = normalise2dpts(data_orig(4:6,:));
data_norm = [data_norm_img1; data_norm_img2];
%outlier removal - Multi-GS (RANSAC).
rng(0);
[ ~,res,~,~] = multigsSampling(100, data_norm, 500, 10);
con = sum(res<=0.05);
[ ~, maxinx ] = max(con);
inliers = find(res(:,maxinx)<=0.05); 

matches_1 = pts1(:, inliers);
matches_2 = pts2(:, inliers);

% delete duplicate feature match
[h1, indh1] = sortrows(matches_1');
[~,  indk1] = unique(h1,'rows');
matches_1 = matches_1(:, indh1(indk1));
matches_2 = matches_2(:, indh1(indk1));
[h2, indh2] = sortrows(matches_2');
[~, indk2] = unique(h2, 'rows');
matches_1 = matches_1(:, indh2(indk2));
matches_2 = matches_2(:, indh2(indk2));

% show the results of RANSAC
% figure;
% imshow([img1 img2]);
% title('Ransac''s results');
% hold on;
% plot(data_orig(1,:),data_orig(2,:),'ro','LineWidth',2);
% plot(data_orig(4,:)+size(img1,2),data_orig(5,:),'ro','LineWidth',2);
% for i=1:length(inliers)
%     plot(data_orig(1,inliers(i)),data_orig(2,inliers(i)),'go','LineWidth',2);
%     plot(data_orig(4,inliers(i))+size(img1,2),data_orig(5,inliers(i)),'go','LineWidth',2);
%     plot([data_orig(1,inliers(i)) data_orig(4,inliers(i))+size(img1,2)],[data_orig(2,inliers(i)) data_orig(5,inliers(i))],'g-');
% end

end