function [ local_sample_lines, local_slope_lines, global_sample_lines, global_slope_lines ] = alllinesDetect( imgpath, img, C1, C2 )
% detect line segments in image using LSD: "A fast line segment detector with a
% false detection control" by R. G. von Gioi, J. Jakubowicz, J.-M. Morel, and G. Randall. Lsd
%% get the start_points and end_points of each straight line use LSD
lines = lsd(imgpath);
[M, N, ~] = size(img);
x_dis = (N-1)/C2;
y_dis = (M-1)/C1;  % width and height of each grid
dia_line = sqrt(x_dis^2+y_dis^2);

%% threshold for all lines
len_threshold = 0.5*dia_line;  % length threshold of line segments, ignore all segments of length less than threshold
len_lines = sqrt( (lines(2,:)-lines(1,:)).^2 + (lines(4,:)-lines(3,:)).^2 );
lines(6,:) = len_lines;
lines_ind = (len_lines>=len_threshold);
lines = lines(:,lines_ind);
[~, index] = sort(lines(6,:));
lines = lines(:, index);

%% plot all the LSD lines.
% figure,imshow(img);
% title(strcat('lines:',32 ,num2str(size(lines,2))));
% hold on
% for k = 1:size(lines, 2)
%     plot(lines(1:2, k), lines(3:4, k), 'LineWidth', 2, 'Color', [1, 0, 0]);
%     text((lines(1,k)+lines(2,k))/2,(lines(3,k)+lines(4,k))/2, num2str(k),'FontSize',8);
% end
% hold off

%% merge the local lines into global lines
mergedlines = mergelinesegments(lines);
if size(mergedlines,2)>0
    mergedlines = mergedlines(:,mergedlines(6,:)>=len_threshold);
end
lines_ = [lines, mergedlines];

%% plot the lines.
% figure,imshow(img);
% title(strcat('global lines and merged lines:',32 ,num2str(size(lines_,2))));
% hold on
% for k = 1:size(lines_, 2)
%     plot(lines_(1:2, k), lines_(3:4, k), 'LineWidth', 2, 'Color', [1, 0, 0]);
%     text((lines_(1,k)+lines_(2,k))/2,(lines_(3,k)+lines_(4,k))/2, num2str(k),'FontSize',8);
% end
% hold off

%% threshold for global lines and local lines
locallines_ind = (lines_(6,:)<2*dia_line);
globallines_ind = (lines_(6,:)>=2*dia_line);

locallines_ = lines_(:,locallines_ind);
globallines_ = lines_(:,globallines_ind);

%% uniformly sample points on lines such that each mesh grid contains a sample point of the line
% local lines
num_lines = max(3, 2*round(locallines_(6,:)./min(x_dis, y_dis)));  % sample numbers of each line segment
local_sample_lines = zeros( 2*size(locallines_,2), max(num_lines)+1);
local_slope_lines = zeros(2*size(locallines_,2),1);  % slope of each lines
for k=1:size(locallines_,2)
    x1 = locallines_(1,k);  y1 = locallines_(3,k);
    x2 = locallines_(2,k);  y2 = locallines_(4,k);
    if x1~=x2
        slope = (y2-y1)/(x2-x1);
        xseq = linspace(x1,x2, num_lines(k));
        yseq = (xseq-x1).*slope+y1;
        r_xseq = (xseq>=1) & (xseq<=N);
        r_yseq = (yseq>=1) & (yseq<=M);
        r_seq = r_xseq & r_yseq;
        local_sample_lines(2*k-1:2*k, 1:sum(r_seq)) = [xseq(r_seq); yseq(r_seq)];
        local_sample_lines(2*k-1:2*k, end) = sum(r_seq); 
        local_slope_lines(2*k-1:2*k) = slope;
    end
    if x1==x2
        xseq = x1.*ones(1, num_lines(k));
        yseq = linspace(y1,y2, num_lines(k));
        r_xseq = (xseq>=1) & (xseq<=N);
        r_yseq = (yseq>=1) & (yseq<=M);
        r_seq = r_xseq & r_yseq;
        local_sample_lines(2*k-1:2*k, 1:sum(r_seq)) = [xseq(r_seq); yseq(r_seq)];
        local_sample_lines(2*k-1:2*k, end) = sum(r_seq);    
        local_slope_lines(2*k-1:2*k) = inf;
    end
end

% global lines
num_lines = max(3, 2*round(globallines_(6,:)./min(x_dis, y_dis)));  % sample numbers of each line segment
global_sample_lines = zeros( 2*size(globallines_,2), max(num_lines)+1);
global_slope_lines = zeros(2*size(globallines_,2),1);  % slope of each lines
for k=1:size(globallines_,2)
    x1 = globallines_(1,k);  y1 = globallines_(3,k);
    x2 = globallines_(2,k);  y2 = globallines_(4,k);
    if x1~=x2
        slope = (y2-y1)/(x2-x1);
        xseq = linspace(x1,x2, num_lines(k));
        yseq = (xseq-x1).*slope+y1;
        r_xseq = (xseq>=1) & (xseq<=N);
        r_yseq = (yseq>=1) & (yseq<=M);
        r_seq = r_xseq & r_yseq;
        global_sample_lines(2*k-1:2*k, 1:sum(r_seq)) = [xseq(r_seq); yseq(r_seq)];
        global_sample_lines(2*k-1:2*k, end) = sum(r_seq); 
        global_slope_lines(2*k-1:2*k) = slope;
    end
    if x1==x2
        xseq = x1.*ones(1, num_lines(k));
        yseq = linspace(y1,y2, num_lines(k));
        r_xseq = (xseq>=1) & (xseq<=N);
        r_yseq = (yseq>=1) & (yseq<=M);
        r_seq = r_xseq & r_yseq;
        global_sample_lines(2*k-1:2*k, 1:sum(r_seq)) = [xseq(r_seq); yseq(r_seq)];
        global_sample_lines(2*k-1:2*k, end) = sum(r_seq);    
        global_slope_lines(2*k-1:2*k) = inf;
    end
end

end

