%输入参数：两条直线的结构体数组
%输出参数：两条直线的交点
function [X, Y]=Itspoint(line1,line2)

        k1=line1.k;
        k2=line2.k;
        b1=line1.b;
        b2=line2.b;
        X=NaN;
        Y=NaN;
    if k1==k2      %两直线平行，无交点。
        %do nothing
    elseif k1~=Inf && k2~=Inf
        if abs(atan((k2 - k1)/(1+ k1*k2)))> 3.14/8  %当两直线夹角<3.14/8，忽略该交点。
            X=(b2-b1)/(k1-k2);
            Y=k1*X+b1;
        end

    elseif k1==Inf       %当一条直线斜率为无穷，另一直线斜率需满足-2.4142<k<2.4142，方估计交点。
        if k2>2.4142 || k2<- 2.4142
        else
            X=line1.point1(1);
            Y=k2*X+b2;
        end
    elseif k2==Inf
        if k1>2.4142 || k1<- 2.4142
        else
            X=line2.point1(1);
            Y=k1*X+b1;
        end
    end

end

