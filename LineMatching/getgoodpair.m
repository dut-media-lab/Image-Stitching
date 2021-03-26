function [ind1,ind2] = getgoodpair(plines1,lines2,dist)

len1=length(plines1);
len2=length(lines2);
ind1=[];
ind2=[];
for i=1:len1
    for j=1:len2
        if   isclose(plines1(i),lines2(j),dist)
            ind1 = [ind1 plines1(i).ind];
            ind2 = [ind2 j];
        end
    end
end


end

function [isc] = isclose(line1,line2,dist)


if disp2line(line1.point1,line2)> dist || disp2line(line1.point2,line2)> dist || ...
        disp2line(line2.point1,line1)> dist || disp2line(line2.point2,line1)> dist || ...
        (norm((line1.point1 + line1.point2)/2 - (line2.point1 + line2.point2)/2 )> ...
        ((norm(line1.point1 - line1.point2)+norm(line2.point1 - line2.point2))/2))
    isc=false;
    return;
else
    isc = true;
    return;
end

end
function dis=disp2line(point,line)
k=line.k;
b=line.b;
if k~=Inf
    dis=abs(k*point(1)-point(2)+b)/sqrt(k*k+1);
else
    dis = abs(point(1)-line.point1(1));
end
end