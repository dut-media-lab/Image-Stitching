function x = vgg_get_nonhomg(x)
% p = vgg_get_nonhomg(h)
%
% Convert a set of homogeneous points to non-homogeneous form
% Points are stored as column vectors, stacked horizontally, e.g.
%  [x0 x1 x2 ... xn ;
%   y0 y1 y2 ... yn ;
%   w0 w1 w2 ... wn ]

% Modified by TW

if isempty(x)
  x = []; 
  return; 
end

d = size(x,1) - 1;
x = x(1:d,:)./(ones(d,1)*x(end,:));

return