function pointlist=intspoints(lines,imsize)
len=length(lines);
k=1;
pointlist(k).point=[0 0];
pointlist(k).lines=[0 0];
for i = 1:len
    for j =i+1:len
        [a,b]=Itspoint(lines(i),lines(j));
        if (a>0) && (a<imsize(2)) && (b>0) && (b<imsize(1))
            
            if  isneighb([a b],lines(i),lines(j))
                pointlist(k).point=[a b];
                pointlist(k).lines=[i j];
                k=k+1;
            end
            
        end
        
    end
end
end

function isneb=isneighb(intp,line1,line2)
isneb=false;
if isnan(intp(1)) || isnan(intp(2))
    return;
end
l1=norm(line1.point1-line1.point2);
l2=norm(line2.point1-line2.point2);
if (norm(intp-line1.point1)<0.2*l1 ||  norm(intp-line1.point2)<0.2*l1 )|| ...
        (norm(intp-line2.point1)<0.2*l2 ||  norm(intp-line2.point2)<0.2*l2)
    isneb=true;
end

end
