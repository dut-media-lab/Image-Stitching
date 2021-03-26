function [ sparse_shortline, sparse_longline ] = energyLineSegments( img, sa_shortlines, sl_shortlines, sa_longlines, sl_longlines, init_H, C1, C2 )
% given detected line segments, devided into short lines, set a strict threshold, and for long lines,set a loose threshold.
%then calculate a sparse matrix to preserving line struture in image separately 
num_V = (C1+1)*(C2+1);  % number of control vertices

X_col = linspace(1,size(img,2),C2+1); % column index of cells
Y_row = linspace(1,size(img,1),C1+1); % row index of cells
x_dis = X_col(2)-X_col(1);  % the width of scale-cell
y_dis = Y_row(2)-Y_row(1);  % the height of scale-cell

Mesh_ps = zeros(4,2); Mesh_pe = zeros(4,2);
%% construct a sparse matrix for short lines
if isempty(sa_shortlines)
    sparse_shortline=[];
else
    row_sps = sum(sa_shortlines(1:2:end-1,end)-1);
    sp_i = zeros(16*row_sps,1); % row index
    sp_j = zeros(16*row_sps,1); % column index
    sp_s = zeros(16*row_sps,1); % value index
    k=1;
    % rotated vertical line equidistant-preserving
    for i=1:2:size(sa_shortlines,1)-1
        num_s = sa_shortlines(i,end);  % number of sample points in this segment
        k_xy = calcSlope(init_H, sl_shortlines(i), [sa_shortlines(i,1); sa_shortlines(i+1,1)]); % line's slope after transformation
        if isinf(abs(k_xy));  nor_vec=[1,0]; else; nor_vec = [k_xy, -1];  end
        nor_vec = nor_vec./norm(nor_vec);  % normal vector of warped lines
        for j=1:num_s-1
            lps = [sa_shortlines(i,j),     sa_shortlines(i+1,j)];
            lpe = [sa_shortlines(i,j+1),   sa_shortlines(i+1,j+1)];
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
            sp_s(16*k-15:16*k) = [-nor_vec(1).*coeff_mesh_ps; nor_vec(1).*coeff_mesh_pe; 
                                  -nor_vec(2).*coeff_mesh_ps; nor_vec(2).*coeff_mesh_pe];
            k = k + 1;
        end
    end

sparse_shortline = sparse(sp_i, sp_j, sp_s, row_sps, 2*num_V);
end
%% construct a sparse matrix for long lines
if isempty(sa_longlines)
    sparse_longline=[];
else
    row_spl = sum(sa_longlines(1:2:end-1,end)-1);
    sp_i = zeros(16*row_spl,1); % row index
    sp_j = zeros(16*row_spl,1); % column index
    sp_s = zeros(16*row_spl,1); % value index
    k=1;
    % rotated vertical line equidistant-preserving
    for i=1:2:size(sa_longlines,1)-1
        num_s = sa_longlines(i,end);  % number of sample points in this segment
        k_xy = calcSlope(init_H, sl_longlines(i), [sa_longlines(i,1); sa_longlines(i+1,1)]); % line's slope after transformation
        if isinf(abs(k_xy));  nor_vec=[1,0]; else; nor_vec = [k_xy, -1];  end
        nor_vec = nor_vec./norm(nor_vec);  % normal vector of warped lines
        for j=1:num_s-1
            lps = [sa_longlines(i,j),     sa_longlines(i+1,j)];
            lpe = [sa_longlines(i,j+1),   sa_longlines(i+1,j+1)];
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
            sp_s(16*k-15:16*k) = [-nor_vec(1).*coeff_mesh_ps; nor_vec(1).*coeff_mesh_pe; 
                                  -nor_vec(2).*coeff_mesh_ps; nor_vec(2).*coeff_mesh_pe];
            k = k + 1;
        end
    end


    sparse_longline = sparse(sp_i, sp_j, sp_s, row_spl, 2*num_V);
end
end