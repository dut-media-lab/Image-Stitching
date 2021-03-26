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


function [pMap, dMap] = doMBS(I, param)

if ~exist('param','var')
    param = getParam();
end

scale = param.MAX_DIM/max(size(I,1),size(I,2));
Ir = imresize(I, scale);

% check the type of the Ir
if ~isa(Ir, 'uint8')
    Ir = mat2gray(Ir);
    Ir = uint8(Ir*255);
end
if size(Ir, 3) == 1
    Ir = repmat(Ir,[1,1,3]);
end

% compute saliency
tic
dMap = MBS(Ir, param.use_lab, param.remove_border, param.use_geodesic);
t1 = toc;


if param.use_backgroundness
    tic
    bgMap = BG(Ir, param.COV_REG,param.MARGIN_RATIO);
    t2 = toc;
    dMap = dMap + bgMap;
end
pMap = dMap;

% postprocess
tic
if param.center_bias
    cmap = imresize(param.cmap, [size(pMap,1), size(pMap,2)]);
    pMap = pMap.*cmap;
end
pMap = mat2gray(pMap);
radius = floor(param.smooth_alpha*sqrt(mean(pMap(:))));
pMap = morphSmooth(pMap, max(radius, 3));
pMap = enhanceContrast(pMap, param.contrast_b);
t3 = toc;

pMap = imresize(pMap, [size(I,1), size(I,2)]);
pMap = mat2gray(pMap);

if param.verbose
    fprintf('computing MBD map: %f\n', t1);
    if param.use_backgroundness
        fprintf('computing BG map: %f\n', t2);
    end
    fprintf('postprocessing: %f\n', t3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sMap = morphSmooth(I,width)
% opening by reconstruction followed by closing by reconstruction
% see the following material for detailed explanations 
% http://www.mathworks.com/products/demos/image/watershed/ipexwatershed.html
I = uint8(I*255);
se = strel('square',width);
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
sMap = mat2gray(Iobrcbr);

function pMap = enhanceContrast(I, b)

t = 0.5*max(I(:));
v1 = mean(I(I>=t));
v2 = mean(I(I<t));
pMap = 1./(1+exp(-b*(I-0.5*(v1+v2)))); 

