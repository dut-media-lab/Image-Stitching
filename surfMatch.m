function [matchpts1,matchpts2]=surfMatch(img1,img2)    
im_ch = size(img1,3);
if im_ch > 1
    gray1 = im2double(rgb2gray(img1));
    gray2 = im2double(rgb2gray(img2));
elseif im_ch == 1
    gray1 = im2double(im1);
    gray2 = im2double(im2);
end
metric_threshold = 200; % use surf feature
num_octaves = 3;
num_scale_levels = 4;
%  ROI = [1 1 size(I,2) size(I,1)]; % Rectangular region of interest,分别表示（x，y，width，height），x，y为区域上角点。
points1 = detectSURFFeatures(gray1, 'MetricThreshold', metric_threshold, 'NumOctaves', num_octaves, 'NumScaleLevels', num_scale_levels);
points2 = detectSURFFeatures(gray2, 'MetricThreshold', metric_threshold, 'NumOctaves', num_octaves, 'NumScaleLevels', num_scale_levels);
    
[features1, valid_points1] = extractFeatures(gray1, points1);
[features2, valid_points2] = extractFeatures(gray2, points2);
    
[indexPairs,matchmetric] = matchFeatures(features1, features2);
    
matched_points1 = valid_points1(indexPairs(:, 1), :);
matched_points2 = valid_points2(indexPairs(:, 2), :);
    
matchpts1 = double(matched_points1.Location');
matchpts2 = double(matched_points2.Location');
% matchpts1 = [X1;ones(1,size(X1, 2))];
% matchpts2 = [X2;ones(1,size(X2, 2))];