function [ptsnear1,ptsnear2]=getnearbypoints(line1,line2,side,img)
[Y,X,~] = size(img);
if side==1
    [ptsnear1, ptsnear2]=getCNpoints(line1,line2,'L',X,Y);
else
    [ptsnear1, ptsnear2]=getCNpoints(line1,line2,'R',X,Y);
end
end

function [pointsnearby1, pointsnearby2] = getCNpoints(line1,line2,side,Xsize,Ysize)
pointsnearby1 = [];
pointsnearby2 = [];
maxp = 40;
d2line = 2;
d2midline = 0.5;

% find the perpendicular bisector
line1_mid = (line1.point1+line1.point2)/2;
line2_mid = (line2.point1+line2.point2)/2;
line1_len = norm(line1.point1-line1.point2);
line2_len = norm(line2.point1-line2.point2);

if line1.k ~= 0
    if line1.k == Inf
        midline1.k = 0;
        midline1.b = line1_mid(2);
    else
        midline1.k = -1/line1.k;
        midline1.b = line1_mid(2)-midline1.k*line1_mid(1);
    end
else
    midline1.k = Inf;
    midline1.b = Inf;
    midline1.point1 = line1_mid;
end
if line2.k ~= 0
    if line2.k == Inf
        midline2.k = 0;
        midline2.b = line2_mid(2);
    else
        midline2.k = -1/line2.k;
        midline2.b = line2_mid(2)-midline2.k*line2_mid(1);
    end
else
    midline2.k = Inf;
    midline2.b = Inf;
    midline2.point1 = line2_mid;
end

if side=='L'
    sidep1=line1.pleft;
    sidep2=line2.pleft;
    if isempty(sidep1) || isempty(sidep2)
        return;
    end
    [~,ind1,ind2] = intersect(sidep1(:,3),sidep2(:,3));
    n = length(ind1);
    if n<5
        return;
    elseif n>maxp
        n = maxp;
        ind1 = ind1(1:maxp);
        ind2 = ind2(1:maxp);
    end
    
    ps1 = sidep1(ind1,1:2);
    ps2 = sidep2(ind2,1:2);
    
    i = n;
    j = n-1;
    for k=1:n-2
        if k==i
            continue;
        end
        K1 = line1.intsect1;
        K2 = line1.intsect2;
        P1 = ps1(i,:);
        P2 = ps1(j,:);
        P3 = ps1(k,:);
        K1P1 = [P1(2)-K1(2), K1(1)-P1(1), P1(1)*K1(2)-K1(1)*P1(2)];
        K1P2 = [P2(2)-K1(2), K1(1)-P2(1), P2(1)*K1(2)-K1(1)*P2(2)];
        K1P3 = [P3(2)-K1(2), K1(1)-P3(1), P3(1)*K1(2)-K1(1)*P3(2)];
        K2P1 = [P1(2)-K2(2), K2(1)-P1(1), P1(1)*K2(2)-K2(1)*P1(2)];
        K2P2 = [P2(2)-K2(2), K2(1)-P2(1), P2(1)*K2(2)-K2(1)*P2(2)];
        K2P3 = [P3(2)-K2(2), K2(1)-P3(1), P3(1)*K2(2)-K2(1)*P3(2)];
        
        K_1 = line2.intsect1;
        K_2 = line2.intsect2;
        P_1 = ps2(i,:);
        P_2 = ps2(j,:);
        P_3 = ps2(k,:);
        K_1P_1 = [P_1(2)-K_1(2), K_1(1)-P_1(1), P_1(1)*K_1(2)-K_1(1)*P_1(2)];
        K_1P_2 = [P_2(2)-K_1(2), K_1(1)-P_2(1), P_2(1)*K_1(2)-K_1(1)*P_2(2)];
        K_1P_3 = [P_3(2)-K_1(2), K_1(1)-P_3(1), P_3(1)*K_1(2)-K_1(1)*P_3(2)];
        K_2P_1 = [P_1(2)-K_2(2), K_2(1)-P_1(1), P_1(1)*K_2(2)-K_2(1)*P_1(2)];
        K_2P_2 = [P_2(2)-K_2(2), K_2(1)-P_2(1), P_2(1)*K_2(2)-K_2(1)*P_2(2)];
        K_2P_3 = [P_3(2)-K_2(2), K_2(1)-P_3(1), P_3(1)*K_2(2)-K_2(1)*P_3(2)];
        
        key1 = solve_cross(K1P1, K2P2);
        key_1 = solve_cross(K_1P_1, K_2P_2);
        if disp2line(key1,line1) < d2line * line1_len && disp2line(key1,midline1) < d2midline * line1_len && sameside(line1,key1,line1_mid+line1.gradient) && ...
                disp2line(key_1,line2) < d2line * line2_len && disp2line(key_1,midline2) < d2midline * line2_len && sameside(line2,key_1,line2_mid+line2.gradient) && ...
                key1(1)<Xsize && key1(1)>0 && key1(2)<Ysize && key1(2)>0 && key_1(1)<Xsize && key_1(1)>0 && key_1(2)<Ysize && key_1(2)>0
            pointsnearby1 = [pointsnearby1,key1'];
            pointsnearby2 = [pointsnearby2,key_1'];
        end
        
        key2 = solve_cross(K1P1, K2P3);
        key_2 = solve_cross(K_1P_1, K_2P_3);
        if disp2line(key2,line1) < d2line * line1_len && disp2line(key2,midline1) < d2midline * line1_len && sameside(line1,key2,line1_mid+line1.gradient) && ...
                disp2line(key_2,line2) < d2line * line2_len && disp2line(key_2,midline2) < d2midline * line2_len && sameside(line2,key_2,line2_mid+line2.gradient) && ...
                key2(1)<Xsize && key2(1)>0 && key2(2)<Ysize && key2(2)>0 && key_2(1)<Xsize && key_2(1)>0 && key_2(2)<Ysize && key_2(2)>0
            pointsnearby1 = [pointsnearby1,key2'];
            pointsnearby2 = [pointsnearby2,key_2'];
        end
        
        key3 = solve_cross(K1P2, K2P1);
        key_3 = solve_cross(K_1P_2, K_2P_1);
        if disp2line(key3,line1) < d2line * line1_len && disp2line(key3,midline1) < d2midline * line1_len && sameside(line1,key3,line1_mid+line1.gradient) && ...
                disp2line(key_3,line2) < d2line * line2_len && disp2line(key_3,midline2) < d2midline * line2_len && sameside(line2,key_3,line2_mid+line2.gradient) && ...
                key3(1)<Xsize && key3(1)>0 && key3(2)<Ysize && key3(2)>0 && key_3(1)<Xsize && key_3(1)>0 && key_3(2)<Ysize && key_3(2)>0
            pointsnearby1 = [pointsnearby1,key3'];
            pointsnearby2 = [pointsnearby2,key_3'];
        end
        
        key4 = solve_cross(K1P2, K2P3);
        key_4 = solve_cross(K_1P_2, K_2P_3);
        if disp2line(key4,line1) < d2line * line1_len && disp2line(key4,midline1) < d2midline * line1_len && sameside(line1,key4,line1_mid+line1.gradient) && ...
                disp2line(key_4,line2) < d2line * line2_len && disp2line(key_4,midline2) < d2midline * line2_len && sameside(line2,key_4,line2_mid+line2.gradient) && ...
                key4(1)<Xsize && key4(1)>0 && key4(2)<Ysize && key4(2)>0 && key_4(1)<Xsize && key_4(1)>0 && key_4(2)<Ysize && key_4(2)>0
            pointsnearby1 = [pointsnearby1,key4'];
            pointsnearby2 = [pointsnearby2,key_4'];
        end
        
        key5 = solve_cross(K1P3, K2P1);
        key_5 = solve_cross(K_1P_3, K_2P_1);
        if disp2line(key5,line1) < d2line * line1_len && disp2line(key5,midline1) < d2midline * line1_len && sameside(line1,key5,line1_mid+line1.gradient) && ...
                disp2line(key_5,line2) < d2line * line2_len && disp2line(key_5,midline2) < d2midline * line2_len && sameside(line2,key_5,line2_mid+line2.gradient) && ...
                key5(1)<Xsize && key5(1)>0 && key5(2)<Ysize && key5(2)>0 && key_5(1)<Xsize && key_5(1)>0 && key_5(2)<Ysize && key_5(2)>0
            pointsnearby1 = [pointsnearby1,key5'];
            pointsnearby2 = [pointsnearby2,key_5'];
        end
        
        key6 = solve_cross(K1P3, K2P2);
        key_6 = solve_cross(K_1P_3, K_2P_2);
        if disp2line(key6,line1) < d2line * line1_len && disp2line(key6,midline1) < d2midline * line1_len && sameside(line1,key6,line1_mid+line1.gradient) && ...
                disp2line(key_6,line2) < d2line * line2_len && disp2line(key_6,midline2) < d2midline * line2_len && sameside(line2,key_6,line2_mid+line2.gradient) && ...
                key6(1)<Xsize && key6(1)>0 && key6(2)<Ysize && key6(2)>0 && key_6(1)<Xsize && key_6(1)>0 && key_6(2)<Ysize && key_6(2)>0
            pointsnearby1 = [pointsnearby1,key6'];
            pointsnearby2 = [pointsnearby2,key_6'];
        end
    end
    
else
    sidep1 = line1.pright;
    sidep2 = line2.pright;
    if isempty(sidep1) || isempty(sidep2)
        return;
    end
    [~,ind1,ind2] = intersect(sidep1(:,3),sidep2(:,3));
    n = length(ind1);
    if n<5
        return;
    elseif n>maxp
        n = maxp;
        ind1 = ind1(1:maxp);
        ind2 = ind2(1:maxp);
    end
    
    ps1 = sidep1(ind1,1:2);
    ps2 = sidep2(ind2,1:2);
    
    i = n;
    j = n-1;
    for k=1:n-2
        if k==i
            continue;
        end
        K1 = line1.intsect1;
        K2 = line1.intsect2;
        P1 = ps1(i,:);
        P2 = ps1(j,:);
        P3 = ps1(k,:);
        K1P1 = [P1(2)-K1(2), K1(1)-P1(1), P1(1)*K1(2)-K1(1)*P1(2)];
        K1P2 = [P2(2)-K1(2), K1(1)-P2(1), P2(1)*K1(2)-K1(1)*P2(2)];
        K1P3 = [P3(2)-K1(2), K1(1)-P3(1), P3(1)*K1(2)-K1(1)*P3(2)];
        K2P1 = [P1(2)-K2(2), K2(1)-P1(1), P1(1)*K2(2)-K2(1)*P1(2)];
        K2P2 = [P2(2)-K2(2), K2(1)-P2(1), P2(1)*K2(2)-K2(1)*P2(2)];
        K2P3 = [P3(2)-K2(2), K2(1)-P3(1), P3(1)*K2(2)-K2(1)*P3(2)];
        
        K_1 = line2.intsect1;
        K_2 = line2.intsect2;
        P_1 = ps2(i,:);
        P_2 = ps2(j,:);
        P_3 = ps2(k,:);
        K_1P_1 = [P_1(2)-K_1(2), K_1(1)-P_1(1), P_1(1)*K_1(2)-K_1(1)*P_1(2)];
        K_1P_2 = [P_2(2)-K_1(2), K_1(1)-P_2(1), P_2(1)*K_1(2)-K_1(1)*P_2(2)];
        K_1P_3 = [P_3(2)-K_1(2), K_1(1)-P_3(1), P_3(1)*K_1(2)-K_1(1)*P_3(2)];
        K_2P_1 = [P_1(2)-K_2(2), K_2(1)-P_1(1), P_1(1)*K_2(2)-K_2(1)*P_1(2)];
        K_2P_2 = [P_2(2)-K_2(2), K_2(1)-P_2(1), P_2(1)*K_2(2)-K_2(1)*P_2(2)];
        K_2P_3 = [P_3(2)-K_2(2), K_2(1)-P_3(1), P_3(1)*K_2(2)-K_2(1)*P_3(2)];
        
        key1 = solve_cross(K1P1, K2P2);
        key_1 = solve_cross(K_1P_1, K_2P_2);
        if disp2line(key1,line1) < d2line * line1_len && disp2line(key1,midline1) < d2midline * line1_len && sameside(line1,key1,line1_mid+line1.gradient) && ...
                disp2line(key_1,line2) < d2line * line2_len && disp2line(key_1,midline2) < d2midline * line2_len && sameside(line2,key_1,line2_mid+line2.gradient) && ...
                key1(1)<Xsize && key1(1)>0 && key1(2)<Ysize && key1(2)>0 && key_1(1)<Xsize && key_1(1)>0 && key_1(2)<Ysize && key_1(2)>0
            pointsnearby1 = [pointsnearby1,key1'];
            pointsnearby2 = [pointsnearby2,key_1'];
        end
        
        key2 = solve_cross(K1P1, K2P3);
        key_2 = solve_cross(K_1P_1, K_2P_3);
        if disp2line(key2,line1) < d2line * line1_len && disp2line(key2,midline1) < d2midline * line1_len && sameside(line1,key2,line1_mid+line1.gradient) && ...
                disp2line(key_2,line2) < d2line * line2_len && disp2line(key_2,midline2) < d2midline * line2_len && sameside(line2,key_2,line2_mid+line2.gradient) && ...
                key2(1)<Xsize && key2(1)>0 && key2(2)<Ysize && key2(2)>0 && key_2(1)<Xsize && key_2(1)>0 && key_2(2)<Ysize && key_2(2)>0
            pointsnearby1 = [pointsnearby1,key2'];
            pointsnearby2 = [pointsnearby2,key_2'];
        end
        
        key3 = solve_cross(K1P2, K2P1);
        key_3 = solve_cross(K_1P_2, K_2P_1);
        if disp2line(key3,line1) < d2line * line1_len && disp2line(key3,midline1) < d2midline * line1_len && sameside(line1,key3,line1_mid+line1.gradient) && ...
                disp2line(key_3,line2) < d2line * line2_len && disp2line(key_3,midline2) < d2midline * line2_len && sameside(line2,key_3,line2_mid+line2.gradient) && ...
                key3(1)<Xsize && key3(1)>0 && key3(2)<Ysize && key3(2)>0 && key_3(1)<Xsize && key_3(1)>0 && key_3(2)<Ysize && key_3(2)>0
            pointsnearby1 = [pointsnearby1,key3'];
            pointsnearby2 = [pointsnearby2,key_3'];
        end
        
        key4 = solve_cross(K1P2, K2P3);
        key_4 = solve_cross(K_1P_2, K_2P_3);
        if disp2line(key4,line1) < d2line * line1_len && disp2line(key4,midline1) < d2midline * line1_len && sameside(line1,key4,line1_mid+line1.gradient) && ...
                disp2line(key_4,line2) < d2line * line2_len && disp2line(key_4,midline2) < d2midline * line2_len && sameside(line2,key_4,line2_mid+line2.gradient) && ...
                key4(1)<Xsize && key4(1)>0 && key4(2)<Ysize && key4(2)>0 && key_4(1)<Xsize && key_4(1)>0 && key_4(2)<Ysize && key_4(2)>0
            pointsnearby1 = [pointsnearby1,key4'];
            pointsnearby2 = [pointsnearby2,key_4'];
        end
        
        key5 = solve_cross(K1P3, K2P1);
        key_5 = solve_cross(K_1P_3, K_2P_1);
        if disp2line(key5,line1) < d2line * line1_len && disp2line(key5,midline1) < d2midline * line1_len && sameside(line1,key5,line1_mid+line1.gradient) && ...
                disp2line(key_5,line2) < d2line * line2_len && disp2line(key_5,midline2) < d2midline * line2_len && sameside(line2,key_5,line2_mid+line2.gradient) && ...
                key5(1)<Xsize && key5(1)>0 && key5(2)<Ysize && key5(2)>0 && key_5(1)<Xsize && key_5(1)>0 && key_5(2)<Ysize && key_5(2)>0
            pointsnearby1 = [pointsnearby1,key5'];
            pointsnearby2 = [pointsnearby2,key_5'];
        end
        
        key6 = solve_cross(K1P3, K2P2);
        key_6 = solve_cross(K_1P_3, K_2P_2);
        if disp2line(key6,line1) < d2line * line1_len && disp2line(key6,midline1) < d2midline * line1_len && sameside(line1,key6,line1_mid+line1.gradient) && ...
                disp2line(key_6,line2) < d2line * line2_len && disp2line(key_6,midline2) < d2midline * line2_len && sameside(line2,key_6,line2_mid+line2.gradient) && ...
                key6(1)<Xsize && key6(1)>0 && key6(2)<Ysize && key6(2)>0 && key_6(1)<Xsize && key_6(1)>0 && key_6(2)<Ysize && key_6(2)>0
            pointsnearby1 = [pointsnearby1,key6'];
            pointsnearby2 = [pointsnearby2,key_6'];
        end
    end
end

end

function [point] = solve_cross(line1,line2)
A1 = line1(1);
B1 = line1(2);
C1 = line1(3);
A2 = line2(1);
B2 = line2(2);
C2 = line2(3);
m = A1*B2-A2*B1;
if m==0
    point = [];
else
    point = [(C2*B1-C1*B2)/m, (C1*A2-C2*A1)/m];
end

end

function dis=disp2line(point,line)
k = line.k;
b = line.b;
if k~=Inf
    dis = abs(k*point(1)-point(2)+b)/sqrt(k*k+1);
else
    dis = abs(point(1)-line.point1(1));
end
end