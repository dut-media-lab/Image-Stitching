function [ lines ] = linegradient( I, lines )

I = double(I);

for i=1:length(lines)
    ldx=0;ldy=0;
    p1=lines(i).point1;p2=lines(i).point2;
    
    k1=lines(i).k;
    b1=lines(i).b;
    if k1>-1 && k1<1
        for ii=min(p1(1),p2(1)):max(p1(1),p2(1))
            jj = round(k1 * ii+b1);
            [dx,dy]=pointgrad(I, ii, jj);
            ldx=ldx+dx;
            ldy=ldy+dy;
        end
    else
        for jj=min(p1(2),p2(2)):max(p1(2),p2(2))
            if k1~=Inf
                ii =  round((jj-b1)/k1);
            else
                ii = p1(1);
            end
            [dx,dy]=pointgrad(I, ii, jj);
            ldx=ldx+dx;
            ldy=ldy+dy;
        end
    end
    
    lines(i).gradient=[ldx,ldy]./max(abs(ldx),abs(ldy))*15;%,atan(ldy/ldx)
end

end
function [dx,dy] = pointgrad(I, i, j)

i=floor(i);
j=floor(j);
s=size(I);
if j<2||j>=s(1)||i<2||i>=s(2)
    dx=0;
    dy=0;
else
    dx = (I(j,i+1) - I(j,i-1))/2;
    dy = (I(j+1,i) - I(j-1,i))/2;
end
end

