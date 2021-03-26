function H = calcHomo(pts1,pts2)

%  H*pts1 = pts2
% Normalise point distribution.
data_pts = [ pts1; ones(1,size(pts1,2)) ; pts2; ones(1,size(pts2,2)) ];
[ dat_norm_pts1,T1 ] = normalise2dpts(data_pts(1:3,:));
[ dat_norm_pts2,T2 ] = normalise2dpts(data_pts(4:6,:));
data_norm = [ dat_norm_pts1 ; dat_norm_pts2 ];

%-----------------------
% Global homography (H).
%-----------------------
%fprintf('DLT (projective transform) on inliers\n');
% Refine homography using DLT on inliers.
%fprintf('> Refining homography (H) using DLT...');
[ h,~,~,~ ] = feval('homography_fit',data_norm);
H = T2\(reshape(h,3,3)*T1);


end