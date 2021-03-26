function [canlines] = addpointsnearby(lines,pointlist,sublinds,charap)

canlines=lines(sublinds);
len=length(sublinds);
for i=1:len
    
    [canlines(i).intsect1, canlines(i).intsect2] = get2imps(sublinds(i),lines,pointlist);
    
    [canlines(i).pleft, canlines(i).pright]=addcharapsrect(lines(sublinds(i)),charap);
end

end

function [p1,p2]=get2imps(lind,lines,pointlist)

subline=lines(lind);
lnum=length(subline);
p1=zeros(lnum,2);
p2=zeros(lnum,2);
[pinds]=getpoints(lind,pointlist);
n=length(pinds);


points=zeros(n,2);
dist1=zeros(1,n);
dist2=zeros(1,n);
for j=1:n
    p=pointlist(pinds(j)).point;
    points(j,:)=p;
    dist1(j)=norm(p-subline.point1);
    dist2(j)=norm(p-subline.point2);
end

[d1,id1]=min(dist1);
[d2,id2]=min(dist2);
if isempty(d1)
    d1=Inf;
end
if isempty(d2)
    d2=Inf;
end
linelen=norm(subline.point1-subline.point2)/10;
if d1>linelen &&d2<=linelen
    p1 = subline.point1;
    p2=points(id2,:);
elseif d2>linelen && d1<=linelen
    p2 = subline.point2;
    p1=points(id1,:);
elseif d1<=linelen &&d2<=linelen
    p1=points(id1,:);
    p2=points(id2,:);
else
    p1 = subline.point1;
    p2 = subline.point2;
end
mid=(subline.point1+subline.point2)/2;
endp= mid + subline.gradient;
cp=crossproduct(p1,endp,mid);
if cp>0  %cp<0
    t=p1;
    p1=p2;
    p2=t;
end
end

function  [pleft, pright] = addcharapsrect(line,charap)
d2line=2;
d2midline=0.5;
pleft=[];
pright=[];
mid=(line.point1+line.point2)/2;
pg=mid+line.gradient;
linelen= norm(line.point1-line.point2);
midline=line;
if line.k ~=0
    midline.k = -1/line.k;
    if line.k == Inf
        midline.b=line.point1(1);
    else
        midline.b=(line.k * mid(2) + mid(1))/line.k;
    end
else
    midline.k=Inf;
    midline.b=Inf;
end
pointnum=size(charap,1);

for i=1:pointnum
    p=charap(i,:);
    if disp2line(p,line) < d2line * linelen && disp2line(p,midline) < d2midline * linelen
        if sameside(line,p,pg)
            pright=[ pright; [p i] ];
        else
            pleft=[ pleft; [p i] ];
        end
    end
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
