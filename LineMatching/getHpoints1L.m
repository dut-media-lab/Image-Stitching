function [p1,p2] = getHpoints1L(line1,line2,side)
minnum=10;
p1=[];p2=[];
if side == 1
    [C,ind1,ind2]=intersect( line1.pleft(:,3), line2.pleft(:,3));
    n=length(C);
    if n>=minnum
        p1= (line1.pleft(ind1,1:2)); %line1.intsect1;line1.intsect2;
        p2= (line2.pleft(ind2,1:2));%line2.intsect1;line2.intsect2;
    end
    
elseif side == 2
    [C,ind1,ind2]=intersect( line1.pright(:,3), line2.pright(:,3));
    n=length(C);
    if n>=minnum
        p1= (line1.pright(ind1,1:2));%line1.intsect1;line1.intsect2;
        p2= (line2.pright(ind2,1:2));%line2.intsect1;line2.intsect2;
    end
end


end