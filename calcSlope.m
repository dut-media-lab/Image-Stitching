function [ slope ] = calcSlope( init_H, k, pts )
% calculate slope of original lines of slope k after transformation (lies through pts)
h = init_H';
x = pts(1); y = pts(2);
if ~isinf(abs(k))
    slope = (k*x*h(4)*h(8)-k*x*h(5)*h(7)+k*h(6)*h(8)-y*h(4)*h(8)+y*h(5)*h(7)-k*h(5)+h(6)*h(7)-h(4))/...
            (k*x*h(1)*h(8)-k*x*h(2)*h(7)+k*h(3)*h(8)-y*h(1)*h(8)+y*h(2)*h(7)-k*h(2)+h(3)*h(7)-h(1));
else
    slope =  (x*h(4)*h(8)-x*h(5)*h(7)+h(6)*h(8)-h(5))/(x*h(1)*h(8)-x*h(2)*h(7)+h(3)*h(8)-h(2));
end

end

