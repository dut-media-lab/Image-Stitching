function [linematch1, linematch2, lpts1, lpts2] = twoLineMatch(imgpath1, imgpath2, pts1, pts2, parameters)

img1 = imread(imgpath1);
img2 = imread(imgpath2);

% detect line segments
lines_1 = lsd(imgpath1);             % lsd (x1; x2; y1; y2; w)
lines_1 = refineLine(lines_1, img1); % refineline (x1; x2; y1; y2)
lines_2 = lsd(imgpath2);
lines_2 = refineLine(lines_2, img2);

% length of line segments
len_lines1 = sqrt( (lines_1(2,:)-lines_1(1,:)).^2 + (lines_1(4,:)-lines_1(3,:)).^2 );
len_lines2 = sqrt( (lines_2(2,:)-lines_2(1,:)).^2 + (lines_2(4,:)-lines_2(3,:)).^2 );

len_threshold = parameters.line_threshold;
lines1 = lines_1(1:4,len_lines1>=len_threshold);  % target
lines2 = lines_2(1:4,len_lines2>=len_threshold);  % reference
lines1 = [lines1(1,:)', lines1(3,:)', lines1(2,:)', lines1(4,:)'];%(x1, y1, x2, y2)
lines2 = [lines2(1,:)', lines2(3,:)', lines2(2,:)', lines2(4,:)'];


fprintf('> Reading files and preparing...');tic;
% read the images and lines, output the lines struct of endpoints, k, b, gradient
[lines1, pointlist1] = paras(img1, lines1);
[lines2, pointlist2] = paras(img2, lines2);

len1 = length(lines1);
len2 = length(lines2);
sublinds1 = 1:length(lines1);
sublinds2 = 1:length(lines2);
lines1 = addpointsnearby(lines1, pointlist1, sublinds1, pts1');
lines2 = addpointsnearby(lines2, pointlist2, sublinds2, pts2');

fprintf('done (%fs)\n',toc);

fprintf('> Calculating similarities between line neighborhoods...');tic;
% calculate the similarities between line neighborhood
simL=zeros(len1,len2);
simR=zeros(len1,len2);

for i=1:len1
    t=lines1(i);
    for j=1:len2
        [simL(i,j),simR(i,j)] = distline(t,lines2(j));
    end
end

k=[];
for i=1:len1
    for j=1:len2
        if simL(i,j)>0.95 && (simL(i,j)==max(simL(i,:)) && simL(i,j)==max(simL(:,j)))
            k = [k;[i,j]];
            break;
        end
    end
end
simside1 = ones(1,size(k,1));

for i=1:len1
    for j=1:len2
        if simR(i,j)>0.95 && (simR(i,j)==max(simR(i,:)) && simR(i,j)==max(simR(:,j)))
            k = [k;[i,j]];
            break;
        end
    end
end
simeside1 = [simside1 2*ones(1,size(k,1))];

len = size(k,1);
votecan = zeros(len2,len1);
fprintf('done (%fs)\n',toc);

fprintf('> Matching lines and getting matchpoints near matchlines ...');tic;
% matching lines and getting matchpoints near matchlines
for i=1:len
    
    [p1, p2]=getHpoints1L(lines1(k(i,1)),lines2(k(i,2)),simeside1(i));
    
    if ~isempty(p1)
        [F1, ~, ~] = estimateGeometricTransform(p1, p2, 'projective');
        plines = projline(F1.T, lines1);
        [ind11, ind12] = getgoodpair(plines, lines2, 3);
        
        plines = projline(inv(F1.T), lines2);
        [ind22, ind21] = getgoodpair(plines, lines1, 3);
        
        if isempty(ind11)||isempty(ind22)
            continue;
        end
        
        [indfinal]=intersect([ind11;ind12]',[ind21;ind22]','rows');
        
        if ~isempty(indfinal)
            indfinal = indfinal';
            ind1 = indfinal(1,:);
            ind2 = indfinal(2,:);
        else
            ind1 = [];ind2=[];
        end
        
        if simeside1(i)==1
            v = simL(k(i,1),k(i,2));
        elseif  simeside1(i)==2
            v = simR(k(i,1),k(i,2));
        end
        
        votecan(ind2+(ind1-1)*len2) = votecan(ind2+(ind1-1)*len2)+v;
    end
end

if size(votecan,1)==1
    [num,ind] = sort(votecan,'descend');
    num = num(1,:);
    ind = ind(1,:);
    votecan = votecan';
    [num2,ind2] = sort(votecan,'descend');
    num2 = num2(1,:);
    ind2 = ind2(1,:);
    k = [];
    if num(1) > 0.9
        k = [k,1];
    end
    nummatch = length(k);
elseif size(votecan,2)==1
    [num,ind] = sort(votecan,'descend');
    num = num(1,:);
    ind = ind(1,:);
    votecan = votecan';
    [num2,ind2] = sort(votecan,'descend');
    num2 = num2(1,:);
    ind2 = ind2(1,:);
    k = [];
    if num2(1) > 0.9
        k = [k,1];
    end
else
    [num,ind] = sort(votecan,'descend');
    num = num(1,:);
    ind = ind(1,:);
    votecan = votecan';
    [num2,ind2] = sort(votecan,'descend');
    num2 = num2(1,:);
    ind2 = ind2(1,:);
    k = [];
    for i=1:length(ind)
        if i==ind2(ind(i)) && (num(i) > 0.9 && num2(ind(i))> 0.9)
            k = [k,i];
        end
    end
    nummatch = length(k);
end

% draw(img1,lines1(k),1:nummatch,num2str(nummatch));
% draw(img2,lines2(ind(k)),1:nummatch,num2str(nummatch));
linestruct1 = lines1(k);
linestruct2 = lines2(ind(k));
linematch1 = zeros(length(linestruct1),4);
linematch2 = linematch1;

ptsnearbyline1=[];
ptsnearbyline2=[];
for i = 1:length(k)
    [lpts1, lpts2] = getnearbypoints(lines1(k(i)), lines2(ind(k(i))),1, img1);
    ptsnearbyline1 = [ptsnearbyline1,lpts1];
    ptsnearbyline2 = [ptsnearbyline2,lpts2];
    [lpts1, lpts2] = getnearbypoints(lines1(k(i)), lines2(ind(k(i))),2, img1);
    ptsnearbyline1 = [ptsnearbyline1,lpts1];
    ptsnearbyline2 = [ptsnearbyline2,lpts2];
end
lpts1 = ptsnearbyline1;
lpts2 = ptsnearbyline2;

for i=1:length(linestruct1)
    linematch1(i,1:2)= linestruct1(i).point1;  % [x1, y1]
    linematch1(i,3:4)= linestruct1(i).point2;  % [x2, y2]
    linematch2(i,1:2)= linestruct2(i).point1;  % [x1, y1]
    linematch2(i,3:4)= linestruct2(i).point2;  % [x2, y2]
end

% delete outliers of line match
linematch11 = [linematch1(:,1),linematch1(:,3), linematch1(:,2), linematch1(:,4)];
linematch22 = [linematch2(:,1),linematch2(:,3), linematch2(:,2), linematch2(:,4)];
linematch1 = refineLine(linematch11', img1);
linematch2 = refineLine(linematch22', img2);

% transform linematch to (N*4), where x1,y1,x2,y2
linematch1([2,3],:) = linematch1([3,2],:);
linematch2([2,3],:) = linematch2([3,2],:);
[linematch1, linematch2] = linesDelete(linematch1', linematch2', pts1, pts2);

% draw(img1,lines1(k),1:nummatch,num2str(nummatch));
% draw(img2,lines2(ind(k)),1:nummatch,num2str(nummatch));
fprintf('done (%fs)\n',toc);

end

% result of linematch
function draw(I, lines, orders, name)
len=length(orders);
figure,
imshow(I),
hold on
if exist('name','var')
    title(strcat('Matching lines: ',32 ,name));
end
for k = 1:len
    if orders(k)~=0
        xy = [lines(orders(k)).point1; lines(orders(k)).point2];
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
        text((xy(1,1)+xy(2,1))/2,(xy(1,2)+xy(2,2))/2,num2str(k));
    end
end
hold off;

end
