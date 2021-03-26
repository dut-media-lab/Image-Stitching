%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Implemetation of the saliency detction method described in paper
%	"Minimum Barrier Salient Object Detection at 80 FPS", Jianming Zhang, 
%	Stan Sclaroff, Zhe Lin, Xiaohui Shen, Brian Price, Radomir Mech, ICCV, 
%       2015
%	
%	Copyright (C) 2015 Jianming Zhang
%
%	This program is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%
%	This program is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%	If you have problems about this software, please contact: 
%       jimmie33@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function param = getParam()

param.MAX_DIM = 300;                        % max image dimension in computation
param.use_lab = true;                       % use Lab color space
param.remove_border = true;                 % detect and remove artificial image frames
param.use_geodesic = false;                 % flag for replacing geodesic distance with MBD
param.use_backgroundness = true;            % MB+
param.COV_REG = 50;                         % covariance regularization term for the maximum value of 255*255
param.MARGIN_RATIO = 0.1;                   % the boundary margion for computing backgroundness map
param.cmap = getCenterMap(param.MAX_DIM);   % the center distance map for center bias
param.center_bias = true;
param.smooth_alpha = 50;                    % see eqn. 9
param.contrast_b = 10;                      % see eqn. 11
param.verbose = false;

function cmap = getCenterMap(dim)
X = [1:dim] - dim/2; X = repmat(X.^2,[dim,1]);
Y = X';
cmap = sqrt(X+Y);
cmap = 1-mat2gray(cmap);