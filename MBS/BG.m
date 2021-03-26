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

function bgMap = BG(I, reg, m_ratio)

I = single(I);
rowMargin = round(m_ratio*size(I,1));
colMargin = round(m_ratio*size(I,2));
mRegion = [];
mRegion{end+1} = I(1:rowMargin,:,:);
mRegion{end+1} = I(end-rowMargin+1:end,:,:);
mRegion{end+1} = I(:,1:colMargin,:);
mRegion{end+1} = I(:,end-colMargin+1:end,:);

tmpMap = zeros(size(I,1)*size(I,2),4);
for i = 1:4
    R = reshape(mRegion{i}, size(mRegion{i},1)*size(mRegion{i},2),[]);
    covMat = cov(R) + reg*eye(size(I,3));
    [U,S,V] = svd(covMat);
    P = U*diag(1./diag(sqrt(S)));
    V = sum(abs(bsxfun(@minus, reshape(I, size(I,1)*size(I,2),[]), mean(R))*P),2);
    tmpMap(:,i) = mat2gray(V);
end

bgMap = sum(tmpMap,2) - max(tmpMap, [], 2);
bgMap = reshape(bgMap,size(I,1),size(I,2));
bgMap = mat2gray(bgMap);

