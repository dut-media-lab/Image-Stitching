function result = gradient_blend(source, mask, target)
% Uses gradient-domain editing to blend a source image into a target image.
% The target image is the background onto which a region of the source
% image is copied and blended; the region of the source image to be
% transferred is given by a binary mask.
%
% Inputs:
%   source: The image from which pixels will be transferred onto the target.
%   mask:   A binary matrix specifying the region of the source image that
%           should be blended into the target.
%   target: The image into which the selected source region is blended;
%           this serves as the background for blending.
% 
% Outputs:
%   result: An image of the same dimensions as the source/target,
%           representing the output of the gradient domain blending.
%
% This function assumes that the inputs, source, mask, and target, all have
% the same width and height.


% To simplify edge cases, we add a 1-pixel border around the source,
% target, and mask. For the source and target images, the extra 1-pixel
% border is created by copying pixel values along the edge, essentially
% extending the original images. For the mask, the extra 1-pixel border
% just consists of all 0s, since we do not want that border to be selected.
% This border is removed after blending.
source = padarray(source, [1,1], 'symmetric');
target = padarray(target, [1,1], 'symmetric');
mask = padarray(mask, [1,1]);

[t_rows, t_cols, ~] = size(target);

% We reshape the source and target to have dimensions t_rows*t_cols x 3, 
% turning each color channel into a column vector. This greatly simplifies
% later computations, which are performed across all color channels
% simultaneously.
s = reshape(source, t_rows*t_cols, []);
t = reshape(target, t_rows*t_cols, []);

% Allocate the RHS vector b.
b = zeros(t_rows*t_cols, 3);

% disp('Constructing the matrix A...');
% tic

% We construct the matrix A efficiently from a set of three vectors:
% row_vec has entries that represent row indexes of A; col_vec has entries
% that represent column indexes of A, and value_vec has entries that
% represent the values at specific positions inside A. These three vectors
% are correlated, such that the final matrix A will have (row_vec(index),
% col_vec(index)) = value_vec(index). Thus, for each entry in A, we add one
% entry into each of row_vec, col_vec, and value_vec.
row_vec = zeros(t_rows*t_cols, 1);
col_vec = zeros(t_rows*t_cols, 1);
value_vec = zeros(t_rows*t_cols, 1);

% Each row of the sparse matrix A represents a linear equation; this
% variable is used to keep track of the current row inside A.
equation_num = 1;

% The matrix A has one equation for each pixel in the target image; we
% iterate through them, and insert the appropriate values in the matrix A
% and the corresponding values in b:
for index = 1:t_rows*t_cols
    if mask(index)
        % 
        b(index,:) = 4*s(index,:) - s(index-1,:) - s(index+1,:) - s(index+t_rows,:) - s(index-t_rows,:);
        
        % Insert a 4 into A at the index of the current central pixel.
        row_vec(equation_num) = index;
        col_vec(equation_num) = index;
        value_vec(equation_num) = 4;
        equation_num = equation_num + 1;
        
        % Insert a -1 for the pixel below the current pixel:
        row_vec(equation_num) = index;
        col_vec(equation_num) = index + 1;
        value_vec(equation_num) = -1;
        equation_num = equation_num + 1;

        % Insert a -1 for the pixel above the current pixel:
        row_vec(equation_num) = index;
        col_vec(equation_num) = index - 1;
        value_vec(equation_num) = -1;
        equation_num = equation_num + 1;

        % Insert a -1 for the pixel to the left of the current pixel:
        row_vec(equation_num) = index;
        col_vec(equation_num) = index - t_rows;
        value_vec(equation_num) = -1;
        equation_num = equation_num + 1;
        
        % Insert a -1 for the pixel to the right of the current pixel:    
        row_vec(equation_num) = index;
        col_vec(equation_num) = index + t_rows;
        value_vec(equation_num) = -1;
        equation_num = equation_num + 1;
    else
        % If the current pixel location is not in the mask, the final value
        % in the blended image is the same as the original value in the
        % target image, so we insert a 1 in the matrix A, and copy the
        % target value into the appropriate position of the vector b:
        row_vec(equation_num) = index;
        col_vec(equation_num) = index;
        value_vec(equation_num) = 1;
        equation_num = equation_num + 1;
        
        b(index,:) = t(index,:);
    end
end

% We create the sparse matrix efficiently:
A = sparse(row_vec, col_vec, value_vec, t_rows*t_cols, t_rows*t_cols);

% toc
% disp('Finished constructing the matrix A...')

% Solve for each color channel:
f_red = A \ b(:,1);
f_green = A \ b(:,2);
f_blue = A \ b(:,3);

% Reshape to the original size:
f_red = reshape(f_red, [t_rows, t_cols]);
f_green = reshape(f_green, [t_rows, t_cols]);
f_blue = reshape(f_blue, [t_rows, t_cols]);

% Stack the channels back together:
result = zeros(t_rows, t_cols, 3);
result(:,:,1) = f_red;
result(:,:,2) = f_green;
result(:,:,3) = f_blue;

% Chop off the border:
result = result(2:t_rows-1, 2:t_cols-1, :);

end