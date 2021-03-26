%% Setup VLFeat toolbox and add other papers' codes.
addpath('modelspecific'); addpath('multigs');  % for feature match and homography

addpath('texture_mapping'); % for our texture mapping

addpath('LSD_matlab'); 

addpath('LineMatching'); % for line segments detection and match

addpath('ASIFT');%for feature point matchlist

addpath('MBS');   % for saliency detection

addpath('LSM');   % for line merge to get long line
addpath('E:\opencv\build\x64\vc12\bin');
% Setup VLFeat toolbox
run('vlfeat-0.9.21/toolbox/vl_setup');