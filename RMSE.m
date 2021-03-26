function [ rmse ] = RMSE( img, C1, C2, pts1, pts2, mesh_X, mesh_Y, off )

X_col = linspace(1,size(img,2), C2+1); % column index of cells
Y_row = linspace(1,size(img,1), C1+1); % row index of cells
x_dis = X_col(2)-X_col(1);  % the width of scale-cell
y_dis = Y_row(2)-Y_row(1);  % the height of scale-cell

Mesh_p = zeros(4, 2);
rmse = 0;

for i=1:size(pts1, 2)
    point = [pts1(1,i), pts1(2,i)];
    px = min( find(point(1)-X_col<x_dis & point(1)-X_col>=0, 1), C2); % the x index of point's position
    py = min( find(point(2)-Y_row<y_dis & point(2)-Y_row>=0, 1), C1); % the y index of point's position
    
    Mesh_p(1:4,:) = [X_col(px), Y_row(py);     % v1
                     X_col(px+1), Y_row(py);   % v2
                     X_col(px+1), Y_row(py+1); % v3
                     X_col(px), Y_row(py+1)];  % v4
    
    coeff_mesh_p = meshGridAlign(Mesh_p, point);
    
    v1 = [mesh_X(py, px), mesh_Y(py, px)];
    v2 = [mesh_X(py, px+1), mesh_Y(py, px+1)];
    v3 = [mesh_X(py+1, px+1), mesh_Y(py+1, px+1)];
    v4 = [mesh_X(py+1, px), mesh_Y(py+1, px)];
    
    warp_point = [coeff_mesh_p(1)*v1(1)+coeff_mesh_p(2)*v2(1)+coeff_mesh_p(3)*v3(1)+coeff_mesh_p(4)*v4(1),
                  coeff_mesh_p(1)*v1(2)+coeff_mesh_p(2)*v2(2)+coeff_mesh_p(3)*v3(2)+coeff_mesh_p(4)*v4(2)];
              
    point_ref = [pts2(1,i); pts2(2,i)];
    rmse = rmse + norm(warp_point - point_ref);
end

rmse = rmse/size(pts1, 2);
rmse = sqrt(rmse);