function [plines] = projline(H,lines)

n=length(lines);
j=1;

for i=1:n
    
    plines(j).point1 = projpoint(H,lines(i).point1);
    plines(j).point2 = projpoint(H,lines(i).point2);
    if (plines(j).point2(1)~=plines(j).point1(1))
        plines(j).k=(plines(j).point2(2)-plines(j).point1(2))/(plines(j).point2(1)-plines(j).point1(1));
        plines(j).b=plines(j).point1(2)-plines(j).k*plines(j).point1(1);
    else
        plines(j).k=Inf;
        plines(j).b=Inf;
    end
    plines(j).ind=i;
    j=j+1;
end

end

function [p]= projpoint(H,point)

p=[point 1];

p=p*H;
p=p./p(3);
p=p(1:2);
end