function [ warped_img ] = meshmap_warp2homo( img, X, Y, wX, wY )
% given original control vertices ori_V and warped control vertices warp_V, 
% map warped image to original image with color information (ÌùÍ¼) use mesh-homography
off = round([ 1 - min(wX(:)) ; 1 - min(wY(:))]);
cw = round(max(wX(:))-min(wX(:))) + 1;
ch = round(max(wY(:))-min(wY(:))) + 1;
C1 = size(X,1)-1;
C2 = size(X,2)-1;

warped_img = texture_mapping_ltl(img, ch, cw, C1, C2, X, Y, wX, wY, off);
clear texture_mapping_ltl;
warped_img = reshape(warped_img,size(warped_img,1),size(warped_img,2)/3,3);
        
end

