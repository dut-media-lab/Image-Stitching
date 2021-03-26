function [ sparse_lineal, c_match ] = energyLineAlign( img, C1, C2, line1, line2 )
%   for given target image, generate C1*C2 mesh grid for line segments align 
%   sparse_lm is a sparse matrix A to solve Ax=0 
% we do not need to delete the all-zero rows in A and crossponding b, as
% A'*A is consistent w or w/o all-zero

[M, N, ~] = size(img);

% line's function: ax+by+c=0
abc_line2 = [line2(:,4)-line2(:,2), line2(:,1)-line2(:,3),...
                line2(:,3).*line2(:,2)-line2(:,1).*line2(:,4)]; 

%%  Generating mesh grid for spatial varying warp
X_col = linspace(1,N,C2+1); % column index of cells
Y_row = linspace(1,M,C1+1); % row index of cells
num_V = (C1+1)*(C2+1);  % number of control vertices 
x_dis = X_col(2)-X_col(1);  % the width of scale-cell
y_dis = Y_row(2)-Y_row(1);  % the height of scale-cell

Mesh_p1 = zeros(4,2); Mesh_p2 = zeros(4,2);
num_line = size(line1,1);% the three indices the sparse function needs
sp_i = ones(16*num_line,1); % row index 
sp_j = ones(16*num_line,1); % column index
sp_s = zeros(16*num_line,1); % value index , ignored zero-element along with i and j
c_match = zeros(2*num_line,1);
%% optimize ||sparse_A*V_star||^2
k=1;
for i=1:num_line
    a=abc_line2(i,1); b=abc_line2(i,2); c=abc_line2(i,3); d=sqrt(a^2+b^2);
    % find the line segments lies in which cell ?
    px1 = min(find( (line1(i,1)-X_col)<x_dis  & (line1(i,1)-X_col)>=0, 1), C2); % the x index of endpoint's position
    py1 = min(find( (line1(i,2)-Y_row)<y_dis  & (line1(i,2)-Y_row)>=0, 1), C1); % the y index of endpoint's position
    px2 = min(find( (line1(i,3)-X_col)<x_dis  & (line1(i,3)-X_col)>=0, 1), C2); % the x index of endpoint's position
    py2 = min(find( (line1(i,4)-Y_row)<y_dis  & (line1(i,4)-Y_row)>=0, 1), C1); % the y index of endpoint's position
    % the cell containing line segments p
    Mesh_p1(1:4,:) = [X_col(px1), Y_row(py1);     % v1
                       X_col(px1+1), Y_row(py1);   % v2
                       X_col(px1+1), Y_row(py1+1); % v3
                       X_col(px1), Y_row(py1+1)];   % v4
    Mesh_p2(1:4,:) = [X_col(px2), Y_row(py2);     % v1
                       X_col(px2+1), Y_row(py2);   % v2
                       X_col(px2+1), Y_row(py2+1); % v3
                       X_col(px2), Y_row(py2+1)];   % v4
    coeff_mesh_p1 = meshGridAlign(Mesh_p1, line1(i,1:2));
    coeff_mesh_p2 = meshGridAlign(Mesh_p2, line1(i,3:4));
    num1 = (C1+1)*(px1-1)+py1; % index of v1*       
    num2 = num1+(C1+1);     % index of v2*        
    num3 = num2+1;            % index of v3*
    num4 = num1+1;            % index of v4*
    num11 = (C1+1)*(px2-1)+py2; % index of v1*       
    num22 = num11+(C1+1);     % index of v2*        
    num33 = num22+1;            % index of v3*
    num44 = num11+1;            % index of v4*    
    sp_i(16*k-15:16*k) = [(2*k-1).*ones(1,8), 2*k.*ones(1,8)];
    sp_j(16*k-15:16*k)= [2*num1-1, 2*num2-1, 2*num3-1, 2*num4-1,...
                       2*num1, 2*num2, 2*num3, 2*num4,...
                       2*num11-1, 2*num22-1, 2*num33-1, 2*num44-1,...
                       2*num11, 2*num22, 2*num33, 2*num44];
    sp_s(16*k-15:16*k)= [a/d.*coeff_mesh_p1; b/d.*coeff_mesh_p1; a/d.*coeff_mesh_p2; b/d.*coeff_mesh_p2];
    c_match(2*k-1:2*k) = [-c/d;-c/d];
    k=k+1;
end

sparse_lineal = sparse(sp_i, sp_j, sp_s, 2*num_line, 2*num_V);

end

