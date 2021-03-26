function [ output_canvas ] = imageBlending( warped_img1, warped_img2, blend_type )
% blend two warped images
w1 = imfill(imbinarize(rgb2gray(warped_img1), 0),'holes');
w2 = imfill(imbinarize(rgb2gray(warped_img2), 0),'holes');

w1 = mat2gray(w1);   w2 = mat2gray(w2);

% composite by average blending
if strcmp(blend_type, 'average')==1
    output_canvas = zeros(size(warped_img1));
    output_canvas(:,:,1) = (warped_img1(:,:,1).*w1 + warped_img2(:,:,1).*w2)./(w1+w2);
    output_canvas(:,:,2) = (warped_img1(:,:,2).*w1 + warped_img2(:,:,2).*w2)./(w1+w2);
    output_canvas(:,:,3) = (warped_img1(:,:,3).*w1 + warped_img2(:,:,3).*w2)./(w1+w2);
end

% composite by linear blending
if strcmp(blend_type, 'linear')==1
    out = warped_img1;
    out_mask = w1;
    % center of out
    [r1, c1] = find(w1);
    out_center1 = [mean(r1) mean(c1)];

    [r2, c2] = find(w2);
    out_center2 = [mean(r2) mean(c2)];
    % compute weighted mask
    vec = out_center2 - out_center1; % vector from out_center to out_i_center
    intsct_mask = w1 & w2; 
    [r, c] = find(intsct_mask(:, :, 1));
    idx = sub2ind(size(w1), r, c);
    out_wmask = zeros(size(w1));
    proj_val = (r - out_center1(1))*vec(1) + (c- out_center1(2))*vec(2); % inner product
    out_wmask(idx) = (proj_val - (min(proj_val)+(1e-3))) / ...
                     ((max(proj_val)-(1e-3)) - (min(proj_val)+(1e-3))); % weight map (of overlapped area) for c1out{i}, 1 channel
    % blending
    mask1 = out_mask(:, :, 1)&(out_wmask==0);
    mask2 = out_wmask;
    mask3 = w2&(out_wmask==0);
    mask1 = cat(3, mask1, mask1, mask1); mask2 = cat(3, mask2, mask2, mask2); mask3 = cat(3, mask3, mask3, mask3);
    output_canvas = out.*(mask1+(1-mask2).*(mask2~=0)) + warped_img2.*(mask2+mask3);
end


end

