function [ rmsedis , rmseerror , rmsecross ] = RMSEline( img, C1, C2, lines1, lines2, mesh_X, mesh_Y, off )
% lines1, lines2 (n*4) matrix, forms (x1,y1,x2,y2)
X_col = linspace(1,size(img,2), C2+1); % column index of cells
Y_row = linspace(1,size(img,1), C1+1); % row index of cells
x_dis = X_col(2)-X_col(1);  % the width of scale-cell
y_dis = Y_row(2)-Y_row(1);  % the height of scale-cell

Mesh_p1 = zeros(4, 2); Mesh_p2 = zeros(4, 2);
rmsedis = 0;
rmseerror = 0;
rmsecross = 0;

for i=1:size(lines1, 1)
    if lines1(i,1)<lines1(i,3)  % 使x轴坐标较小的点为start_point
        start_point = [lines1(i,1), lines1(i,2)];
        end_point = [lines1(i,3), lines1(i,4)];
    else
        end_point = [lines1(i,1), lines1(i,2)];
        start_point = [lines1(i,3), lines1(i,4)];
    end
    
    if start_point(1) ~= end_point(1)
        line_k = (end_point(2)-start_point(2))/(end_point(1)-start_point(1));
        line_b = start_point(2)-line_k*start_point(1);
    else
        line_k = Inf;
        line_b = Inf;
    end
    % sampling
    if line_k ~= Inf
        sample_X = linspace(start_point(1), end_point(1), (end_point(1)-start_point(1))/3+1);
        sample_Y = line_k.*sample_X+line_b;
    else
        if start_point(2)<end_point(2)
            sample_Y = linspace(start_point(2), end_point(2), (end_point(2)-start_point(2))/3+1);
            sample_X = start_point(1).*ones(1, size(sample_Y, 2));
        else
            sample_Y = linspace(end_point(2), start_point(2), (start_point(2)-end_point(2))/3+1);
            sample_X = start_point(1).*ones(1, size(sample_Y, 2));
        end
    end
    sample = [sample_X;sample_Y];  % (2*n)矩阵，n为采样点个数，每列为一个采样点
    
    pxstart = min( find(start_point(1)-X_col<x_dis & start_point(1)-X_col>=0, 1), C2); % the x index of start_point's position
    pystart = min( find(start_point(2)-Y_row<y_dis & start_point(2)-Y_row>=0, 1), C1); % the y index of start_point's position
    pxend = min( find(end_point(1)-X_col<x_dis & end_point(1)-X_col>=0, 1), C2); % the x index of end_point's position
    pyend = min( find(end_point(2)-Y_row<y_dis & end_point(2)-Y_row>=0, 1), C1); % the y index of end_point's position
    
    Mesh_p1(1:4,:) = [X_col(pxstart), Y_row(pystart);     % v1
                     X_col(pxstart+1), Y_row(pystart);   % v2
                     X_col(pxstart+1), Y_row(pystart+1); % v3
                     X_col(pxstart), Y_row(pystart+1)];  % v4
    Mesh_p2(1:4,:) = [X_col(pxend), Y_row(pyend);     % v1
                     X_col(pxend+1), Y_row(pyend);   % v2
                     X_col(pxend+1), Y_row(pyend+1); % v3
                     X_col(pxend), Y_row(pyend+1)];  % v4
    
    coeff_mesh_p1 = meshGridAlign(Mesh_p1, start_point);
    coeff_mesh_p2 = meshGridAlign(Mesh_p2, end_point);
    
    v1start = [mesh_X(pystart, pxstart), mesh_Y(pystart, pxstart)];
    v2start = [mesh_X(pystart, pxstart+1), mesh_Y(pystart, pxstart+1)];
    v3start = [mesh_X(pystart+1, pxstart+1), mesh_Y(pystart+1, pxstart+1)];
    v4start = [mesh_X(pystart+1, pxstart), mesh_Y(pystart+1, pxstart)];
    v1end = [mesh_X(pyend, pxend), mesh_Y(pyend, pxend)];
    v2end = [mesh_X(pyend, pxend+1), mesh_Y(pyend, pxend+1)];
    v3end = [mesh_X(pyend+1, pxend+1), mesh_Y(pyend+1, pxend+1)];
    v4end = [mesh_X(pyend+1, pxend), mesh_Y(pyend+1, pxend)];
    
    warp_start_point = [coeff_mesh_p1(1)*v1start(1)+coeff_mesh_p1(2)*v2start(1)+coeff_mesh_p1(3)*v3start(1)+coeff_mesh_p1(4)*v4start(1),...
                        coeff_mesh_p1(1)*v1start(2)+coeff_mesh_p1(2)*v2start(2)+coeff_mesh_p1(3)*v3start(2)+coeff_mesh_p1(4)*v4start(2)];
    warp_end_point = [coeff_mesh_p2(1)*v1end(1)+coeff_mesh_p2(2)*v2end(1)+coeff_mesh_p2(3)*v3end(1)+coeff_mesh_p2(4)*v4end(1),...
                        coeff_mesh_p2(1)*v1end(2)+coeff_mesh_p2(2)*v2end(2)+coeff_mesh_p2(3)*v3end(2)+coeff_mesh_p2(4)*v4end(2)];                
        
    if lines2(i,1)<lines2(i,3)  % 使x轴坐标较小的点为start_point_ref
        start_point_ref = [lines2(i,1), lines2(i,2)];
        end_point_ref = [lines2(i,3), lines2(i,4)];
    else
        end_point_ref = [lines2(i,1), lines2(i,2)];
        start_point_ref = [lines2(i,3), lines2(i,4)];
    end
    
    line_ref_A = end_point_ref(2)-start_point_ref(2);
    line_ref_B = start_point_ref(1)-end_point_ref(1);
    line_ref_C = end_point_ref(1)*start_point_ref(2)-start_point_ref(1)*end_point_ref(2);
    warp_start_point_dis = abs(line_ref_A*warp_start_point(1)+line_ref_B*warp_start_point(2)+line_ref_C)/sqrt(line_ref_A*line_ref_A+line_ref_B*line_ref_B);
    warp_end_point_dis = abs(line_ref_A*warp_end_point(1)+line_ref_B*warp_end_point(2)+line_ref_C)/sqrt(line_ref_A*line_ref_A+line_ref_B*line_ref_B);
    rmsedis = rmsedis + ((warp_start_point_dis+warp_end_point_dis)/2.0)^2;
    
    line_warp_norm_vector = (warp_end_point-warp_start_point);
    line_ref_norm_vector = (end_point_ref-start_point_ref);
    cp = cross_product(line_warp_norm_vector, line_ref_norm_vector);
    rmsecross = rmsecross+(abs(cp))^2;
    
    Mesh_psample = zeros(4, 2);
    warp_sample = zeros(2, size(sample, 2));
    for j=1:size(sample, 2)
        point = [sample(1,j), sample(2,j)];
        px = min( find(point(1)-X_col<x_dis & point(1)-X_col>=0, 1), C2); % the x index of point's position
        py = min( find(point(2)-Y_row<y_dis & point(2)-Y_row>=0, 1), C1); % the y index of point's position
        
        Mesh_psample(1:4,:) = [X_col(px), Y_row(py);     % v1
                               X_col(px+1), Y_row(py);   % v2
                               X_col(px+1), Y_row(py+1); % v3
                               X_col(px), Y_row(py+1)];  % v4
                           
        coeff_mesh_p = meshGridAlign(Mesh_psample, point);
        
        v1 = [mesh_X(py, px), mesh_Y(py, px)];
        v2 = [mesh_X(py, px+1), mesh_Y(py, px+1)];
        v3 = [mesh_X(py+1, px+1), mesh_Y(py+1, px+1)];
        v4 = [mesh_X(py+1, px), mesh_Y(py+1, px)];
        
        warp_point = [coeff_mesh_p(1)*v1(1)+coeff_mesh_p(2)*v2(1)+coeff_mesh_p(3)*v3(1)+coeff_mesh_p(4)*v4(1),...
                      coeff_mesh_p(1)*v1(2)+coeff_mesh_p(2)*v2(2)+coeff_mesh_p(3)*v3(2)+coeff_mesh_p(4)*v4(2)];
        warp_sample(:, j) = warp_point';                
    end
    
    leastsqu = polyfit(warp_sample(1,:),warp_sample(2,:),1);
    observation = leastsqu(1).*warp_sample(1,:)+leastsqu(2);
    error = abs(observation-warp_sample(2,:));
    sum_error = sum(error)/size(sample, 2);
    rmseerror = rmseerror + sum_error^2;
    
%     line_warp_norm_vector = (warp_end_point-warp_start_point)/norm(warp_end_point-warp_start_point);
%     line_ref_norm_vector = [1,leastsqu(1)]/(1+leastsqu(1)^2);
%     cp = cross_product(line_warp_norm_vector, line_ref_norm_vector);
%     rmsecross = rmsecross+abs(cp);
end

rmsedis = rmsedis/size(lines1, 1);
rmsedis = sqrt(rmsedis);
rmseerror = rmseerror/size(lines1, 1);
rmseerror = sqrt(rmseerror);
rmsecross = rmsecross/size(lines1, 1);
rmsecross = sqrt(rmsecross);
% rmsecross = asind(rmsecross);
end

function cp=cross_product(v1,v2)
x1=v1(1);y1=v1(2);x2=v2(1);y2=v2(2);
cp=x1*y2-x2*y1;
end