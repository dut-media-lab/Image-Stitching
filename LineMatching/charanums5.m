function  CN=charanums5(a,b,c,d,e,f,g,h,i,j)
CN=((a*d - b*c - a*j + b*i + c*j - d*i)*(a*f - b*e - a*h + b*g + e*h - f*g)*(c*h - d*g - c*j + d*i + g*j - h*i))/...
((a*d - b*c - a*h + b*g + c*h - d*g)*(a*h - b*g - a*j + b*i + g*j - h*i)*(c*f - d*e - c*j + d*i + e*j - f*i));
end

% function  CN=charanums5(p5)
% p5=[p5,[1;1;1;1;1]];
% cn51=1-det(p5([1,3,5],:))*det(p5([2,4,5],:))/det(p5([1,4,5],:))/det(p5([2,3,5],:));
% cn52=det(p5([2,3,4],:))*det(p5([1,4,5],:))/det(p5([1,3,4],:))/det(p5([2,4,5],:))-1;
% CN=cn51/cn52;
% end
% 