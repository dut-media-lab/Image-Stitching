function [ refine_lines ] = refineLine( lines, img )
% refine lines in image coordinates,
% lines: 5*M matrix M lines
% |x1,x2,y1,y2|;...
[sz1, sz2, ~] = size(img);
refine_lines = lines;
for i=1:size(lines,2)
    x1=lines(1,i);
    x2=lines(2,i);
    y1=lines(3,i);
    y2=lines(4,i);
    if x1>=1 && x1<=sz2 && x2>=1 && x2<=sz2 && y1>=1 && y1<=sz1 && y2>=1 && y2<=sz1
        continue;        
    end
    a=y2-y1; b=x1-x2; c=x2*y1-x1*y2;
    if abs(b)<=eps  % x1=x2
        x1=min(max(1.1, x1),sz2);
        x2=min(max(1.1, x2),sz2);
        y1=min(max(1.1, y1),sz1);
        y2=min(max(1.1, y2),sz1);
    else
        x1=min(max(1.1, x1),sz2);
        x2=min(max(1.1, x2),sz2);
        y1=-(a*x1+c)/b;
        y2=-(a*x2+c)/b;
    end
    if abs(a)>eps
        y1=min(max(1.1, y1),sz1);
        y2=min(max(1.1, y2),sz1);
        x1=-(b*y1+c)/a;
        x2=-(b*y2+c)/a;
    end
    refine_lines(1,i)=x1; refine_lines(3,i)=y1;
    refine_lines(2,i)=x2; refine_lines(4,i)=y2;
end

end

