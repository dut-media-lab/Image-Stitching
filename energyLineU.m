function [ sparse_us, sparse_ue ] = energyLineU( img, C1, C2, lines_us, lines_ue, init_H )
%  for given target image, C1*C2 mesh grid and lines in target image,
%  generate sparse matrix for line-preserving
% for lines_u, preserve slope-invariant
num_V = (C1+1)*(C2+1);  % number of control vertices
X_col = linspace(1,size(img,2), C2+1); % column index of cells
Y_row = linspace(1,size(img,1), C1+1); % row index of cells
x_dis = X_col(2)-X_col(1);  % the width of scale-cell
y_dis = Y_row(2)-Y_row(1);  % the height of scale-cell

Mesh_ps = zeros(4,2); Mesh_pe = zeros(4,2); Mesh_pm = zeros(4,2);
%% rotated horizontal line, slope-preserving 
hor_sp = sum(lines_us(1:2:end-1,end)-1);
sp_i = zeros(16*hor_sp,1); % row index
sp_j = zeros(16*hor_sp,1); % column index
sp_s = zeros(16*hor_sp,1); % value index
k=1;  init_k = init_H(6)/init_H(3);
for i=1:2:size(lines_us,1)-1
    num_u = lines_us(i,end);
    if num_u<=1; continue; end % if sample points less than 2, continue
    k_xy = calcSlope(init_H, init_k, [lines_us(i,1); lines_us(i+1, 1)]); % slope of u-line after transformation
    if ~isinf(abs(k_xy)); nor_vec = [k_xy, -1]; else; nor_vec = [1, 0]; end
    nor_vec = nor_vec./norm(nor_vec);  % normal vector of warped u-lines
    for j=1:num_u-1
        lps = [lines_us(i,j),     lines_us(i+1,j)];
        lpe = [lines_us(i,j+1),   lines_us(i+1,j+1)];
        pxs = min( find(lps(1)-X_col<x_dis & lps(1)-X_col>=0, 1), C2); % the x index of p's position
        pys = min( find(lps(2)-Y_row<y_dis & lps(2)-Y_row>=0, 1), C1); % the y index of p's position
        pxe = min( find(lpe(1)-X_col<x_dis & lpe(1)-X_col>=0, 1), C2); % the x index of p's position
        pye = min( find(lpe(2)-Y_row<y_dis & lpe(2)-Y_row>=0, 1), C1); % the y index of p's position
        
        nums1 = (C1+1)*(pxs-1) + pys; % index of v1*   
        nums2 = nums1 + C1+1;
        nums3 = nums2 + 1;
        nums4 = nums1 + 1;
        nume1 = (C1+1)*(pxe-1) + pye;
        nume2 = nume1 + C1+1;
        nume3 = nume2 + 1;
        nume4 = nume1 + 1;    
        
        Mesh_ps(1:4,:) = [X_col(pxs), Y_row(pys);     % v1
                            X_col(pxs+1), Y_row(pys);   % v2
                            X_col(pxs+1), Y_row(pys+1); % v3
                            X_col(pxs), Y_row(pys+1)];   % v4
        Mesh_pe(1:4,:) = [X_col(pxe), Y_row(pye);     % v1
                            X_col(pxe+1), Y_row(pye);   % v2
                            X_col(pxe+1), Y_row(pye+1); % v3
                            X_col(pxe), Y_row(pye+1)];   % v4                  
        coeff_mesh_ps = meshGridAlign(Mesh_ps, lps);
        coeff_mesh_pe = meshGridAlign(Mesh_pe, lpe);
        
        sp_i(16*k-15:16*k) = k.*ones(1,16);
        sp_j(16*k-15:16*k) = [2*nums1-1, 2*nums2-1, 2*nums3-1, 2*nums4-1,...
                               2*nume1-1, 2*nume2-1, 2*nume3-1, 2*nume4-1,...
                               2*nums1, 2*nums2, 2*nums3, 2*nums4,...
                               2*nume1, 2*nume2, 2*nume3, 2*nume4];
        sp_s(16*k-15:16*k) = [ -nor_vec(1).*coeff_mesh_ps;  nor_vec(1).*coeff_mesh_pe; -nor_vec(2).*coeff_mesh_ps;  nor_vec(2).*coeff_mesh_pe ];
        k = k + 1;
    end
end
sparse_us = sparse( sp_i, sp_j, sp_s, hor_sp, 2*num_V);

%% rotated horizontal line, equidistant-preserving 
hor_sp = 2*sum(max(0,lines_ue(1:2:end-1,end)-2));
sp_i = zeros(12*hor_sp, 1); % row index
sp_j = zeros(12*hor_sp, 1); % column index
sp_s = zeros(12*hor_sp, 1); % value index
k=1;
for i=1:2:size(lines_ue,1)-1
    num_u = lines_ue(i,end);
    if num_u<=2; continue; end % if sample points less than 2, continue
    for j=2:num_u-1
        lps = [lines_ue(i,j-1),  lines_ue(i+1,j-1)];
        lpm = [lines_ue(i,j),    lines_ue(i+1,j)];
        lpe = [lines_ue(i,j+1),  lines_ue(i+1,j+1)];
        pxs = min( find(lps(1)-X_col<x_dis & lps(1)-X_col>=0, 1), C2); % the x index of p's position
        pys = min( find(lps(2)-Y_row<y_dis & lps(2)-Y_row>=0, 1), C1); % the y index of p's position
        pxm = min( find(lpm(1)-X_col<x_dis & lpm(1)-X_col>=0, 1), C2); % the x index of p's position
        pym = min( find(lpm(2)-Y_row<y_dis & lpm(2)-Y_row>=0, 1), C1); % the y index of p's position
        pxe = min( find(lpe(1)-X_col<x_dis & lpe(1)-X_col>=0, 1), C2); % the x index of p's position
        pye = min( find(lpe(2)-Y_row<y_dis & lpe(2)-Y_row>=0, 1), C1); % the y index of p's position
        
        nums1 = (C1+1)*(pxs-1) + pys; % index of v1*   
        nums2 = nums1 + C1+1;
        nums3 = nums2 + 1;
        nums4 = nums1 + 1;
        numm1 = (C1+1)*(pxm-1) + pym; % index of v1*   
        numm2 = numm1 + C1+1;
        numm3 = numm2 + 1;
        numm4 = numm1 + 1;
        nume1 = (C1+1)*(pxe-1) + pye;
        nume2 = nume1 + C1+1;
        nume3 = nume2 + 1;
        nume4 = nume1 + 1;    
        
        Mesh_ps(1:4,:) = [X_col(pxs), Y_row(pys);     % v1
                            X_col(pxs+1), Y_row(pys);   % v2
                            X_col(pxs+1), Y_row(pys+1); % v3
                            X_col(pxs), Y_row(pys+1)];   % v4
        Mesh_pm(1:4,:) = [X_col(pxm), Y_row(pym);     % v1
                            X_col(pxm+1), Y_row(pym);   % v2
                            X_col(pxm+1), Y_row(pym+1); % v3
                            X_col(pxm), Y_row(pym+1)];   % v4                       
        Mesh_pe(1:4,:) = [X_col(pxe), Y_row(pye);     % v1
                            X_col(pxe+1), Y_row(pye);   % v2
                            X_col(pxe+1), Y_row(pye+1); % v3
                            X_col(pxe), Y_row(pye+1)];   % v4                  
        coeff_mesh_ps = meshGridAlign(Mesh_ps, lps);
        coeff_mesh_pm = meshGridAlign(Mesh_pm, lpm);
        coeff_mesh_pe = meshGridAlign(Mesh_pe, lpe);
        
        sp_i(24*k-23:24*k) = [(2*k-1).*ones(1,12), 2*k.*ones(1,12)];
        sp_j(24*k-23:24*k) = [2*nums1-1, 2*nums2-1, 2*nums3-1, 2*nums4-1,...
                              2*numm1-1, 2*numm2-1, 2*numm3-1, 2*numm4-1,...
                              2*nume1-1, 2*nume2-1, 2*nume3-1, 2*nume4-1,...
                              2*nums1, 2*nums2, 2*nums3, 2*nums4,...
                              2*numm1, 2*numm2, 2*numm3, 2*numm4,...
                              2*nume1, 2*nume2, 2*nume3, 2*nume4];
        sp_s(24*k-23:24*k) = [coeff_mesh_ps; -2.*coeff_mesh_pm; coeff_mesh_pe;...
                                coeff_mesh_ps; -2.*coeff_mesh_pm; coeff_mesh_pe];
        k = k + 1;
    end
end
sparse_ue = sparse( sp_i, sp_j, sp_s, hor_sp, 2*num_V);

end