function pc = z(p,C)
% function pc = vgg_condition_2d(p,C);
%
% Condition a set of 2D homogeneous or nonhomogeneous points using conditioner C

[r,c] = size(p);
if r == 2
  pc = vgg_get_nonhomg(C * vgg_get_homg(p));
elseif r == 3
  pc = C * p;
else
  error ('rows != 2 or 3');
end

end