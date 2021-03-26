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

function [ lines, lengths,angles] = sort_lines_lengths( lines,lengths,angles)
[~,sorted_inds]=sort(lengths,'descend');
lines=lines(sorted_inds,:);
lengths=lengths(sorted_inds);
angles=angles(sorted_inds);
end

