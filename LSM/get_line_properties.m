% This file is part of LSM.
% 
% LSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU v3.0 General Public License as published by
% the Free Software Foundation.
% 
% LSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details
% <http://www.gnu.org/licenses/>.

function [lengths,angles]=get_line_properties(lines)
x1=lines(:,1); x2=lines(:,3);
y1=lines(:,2); y2=lines(:,4);
dx=x2-x1;
dy=y2-y1;
lengths=sqrt(dx.^2+dy.^2);
angles=atan2(dy,dx);
angles(angles<=0)=angles(angles<=0)+pi; %convert from (-pi,pi] to (0,pi]


