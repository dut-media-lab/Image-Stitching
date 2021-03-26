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

function [ lines ] = sort_line_points_y( lines )
% SORT_LINE_POINTS_Y sorts points of lines in a way that points having y
% value greater lies at 1,2 and point having lower y values lies at 3,4
% forcing lower end of each boundary line to be its starting point

for i=1:size(lines)
    tmp_line=lines(i,:);
    if lines(i,2)>lines(i,4)
        %left_line_t=left_line;
    else
        lines(i,:)=[tmp_line(3) tmp_line(4) tmp_line(1) tmp_line(2) tmp_line(5)];
    end    
end

end

