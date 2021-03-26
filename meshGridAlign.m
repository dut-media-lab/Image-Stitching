function [ coeff_Mesh ] = meshGridAlign( Mesh, f_pts )
% 输入四个控制网格点坐标和内部特征点坐标，将特征点表示为四个网格点的线性加权和（双线性插值）
% 返回四个控制网格点线性表示的系数w1,w2,w3,w4  p=w1*v1+w2*v2+w3*v3+w4*v4; 
% f_pts = coeff_Mesh'*Mesh;
%  Mesh: 4*2 matrix (x,y)  v1 ________ v2
%                            |  .p    |
%                            |        |
%                          v4|________| v3
%feature_pts: p = (px,py)
coeff_Mesh = zeros(4,1); % w1 ,w2, w3, w4
area_Mesh = (Mesh(2,1)-Mesh(1,1))*(Mesh(4,2)-Mesh(1,2));  % (v2x-v1x)*(v4y-v1y)

coeff_Mesh(1) = (Mesh(3,1)-f_pts(1))*(Mesh(3,2)-f_pts(2));  % w1 = (v3x-px)*(v3y-py)
coeff_Mesh(2) = (f_pts(1)-Mesh(4,1))*(Mesh(4,2)-f_pts(2));  % w2 = (px-v4x)*(v4y-py)
coeff_Mesh(3) = (f_pts(1)-Mesh(1,1))*(f_pts(2)-Mesh(1,2));  % w3 = (px-v1x)*(py-v1y)
coeff_Mesh(4) = (Mesh(2,1)-f_pts(1))*(f_pts(2)-Mesh(2,2));  % w4 = (v2x-px)*(py-v2y)

coeff_Mesh = coeff_Mesh./area_Mesh;  % coefficient divide area = (v2x-v1x)*(v4y-v1y)

end

