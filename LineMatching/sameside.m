function ss=sameside(line,point1,point2)

ss = false;
if isnan(point2(1)) || isnan(point2(2))
    return;
end

if (line.k~=Inf)
    s1=line.k*point1(1)+line.b-point1(2);
    s2=line.k*point2(1)+line.b-point2(2);
else
    s1=point1(1)-line.point1(1);
    s2=point2(1)-line.point1(1);
end

if s1*s2>0
    ss = true;
end

end
