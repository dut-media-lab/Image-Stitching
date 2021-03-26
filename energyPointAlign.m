function [ sparse_al, psMatch ] = energyPointAlign( img, C1, C2, pts1, pts2 )
%   for given target image, generate C1*C2 mesh grid for feature align 
%   the new alignment term adds scale operator into the sparse matrix
%   sparse_align is a sparse matrix A to solve Ax=b 
% we do not need to delete the all-zero rows  in A and crossponding b, as
% A'*A is equal w or w/o all-zero

scale = 2; % the scale operator for alignment, scale: 1-scale
[M,N,~] = size(img);
%%  Generating mesh grid for spatial varying warp
X_col = linspace(1,N,C2+1); % column index of cells
Y_row = linspace(1,M,C1+1); % row index of cells
num_V = (C1+1)*(C2+1);  % number of control vertices 
x_dis = X_col(2)-X_col(1);  % the width of scale-cell
y_dis = Y_row(2)-Y_row(1);  % the height of scale-cell
max_length = sum((1:scale).^2);

cell_sparse = cell(scale,1);
psMatch = zeros(max_length*2*size(pts1,2), 1);

start = 1;  Mesh_p = zeros(4,2);
for s=1:scale
    % for every feature point, find its position (lie into which cell ?)  
    num_spts = 2*s^2*size(pts1,2); % the three indices the sparse function needs
    sp_i = ones(4*num_spts,1); % row index 
    sp_j = ones(4*num_spts,1); % column index
    sp_s = zeros(4*num_spts,1); % value index , ignored zero-element along with i and j
    pmatch = zeros(num_spts,1);  % feature matches in reference image
    %% optimize ||sparse_A*V_star-Pmatch_2||^2
    k=1;
    for i=1:size(pts1,2)
        % find the feature point lies in which cell ?
        px = min(find( (pts1(1,i)-X_col)<x_dis  & (pts1(1,i)-X_col)>=0, 1), C2); % the x index of p's position
        py = min(find( (pts1(2,i)-Y_row)<y_dis  & (pts1(2,i)-Y_row)>=0, 1), C1); % the y index of p's position
        % count all the scale-cell
        for xi=px-s+1:px
            for yi=py-s+1:py
                if xi>0 && yi>0 && xi+s<=C2+1 && yi+s<=C1+1  % if the scale-cell is well defined  
                    % the cell containing feature p
                    Mesh_p(1:4,:) = [X_col(xi), Y_row(yi);     % v1
                                     X_col(xi+s), Y_row(yi);   % v2
                                     X_col(xi+s), Y_row(yi+s); % v3
                                     X_col(xi), Y_row(yi+s)];   % v4
                    coeff_mesh_p = meshGridAlign(Mesh_p, pts1(:,i));

                    num1 = (C1+1)*(xi-1)+yi; % index of v1*       
                    num2 = num1+s*(C1+1);     % index of v2*        
                    num3 = num2+s;            % index of v3*
                    num4 = num1+s;            % index of v4*
                    sp_i(8*k-7:8*k)=[(2*k-1).*ones(1,4), 2*k.*ones(1,4)];
                    sp_j(8*k-7:8*k)=[2*num1-1, 2*num2-1, 2*num3-1, 2*num4-1,...
                                       2*num1, 2*num2, 2*num3, 2*num4];
                    sp_s(8*k-7:8*k)=[coeff_mesh_p; coeff_mesh_p ];
                    pmatch(2*k-1:2*k) = pts2(1:2,i);
                    k = k+1;
                end
            end
        end
    end
    sp_s = sqrt(1/double(s)).*sp_s; % divide scale operator
    pmatch = sqrt(1/double(s)).*pmatch; 
    sparse_sa = sparse(sp_i, sp_j, sp_s, num_spts, 2*num_V);
    tmp_length = length(pmatch);
    psMatch(start:start-1+tmp_length,:) = pmatch;  % the alignment feature (reference) under scale s
    start = tmp_length + start;
    cell_sparse{s} = sparse_sa;  % the alignment term under scale s
end

sparse_al = [];
for i=1:scale
    sparse_al = [sparse_al; cell_sparse{i}];
end

end

