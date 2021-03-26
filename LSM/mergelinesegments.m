function [mergelines] = mergelinesegments(lines)
% merge linesegments in image coordinates,
% lines: 4*M matrix M lines
% |x1,x2,y1,y2|;...
%step1：对一条直线和另一条直线的两个端点，分别用第一条直线的两个点和第二条直线的第一个点、第二个点做叉积判断是否近似为1.
mergelines=[];
for i=1:size(lines,2)
    for j=i:size(lines,2)
        xi1=lines(1,i); xi2=lines(2,i);
        yi1=lines(3,i); yi2=lines(4,i);
        xj1=lines(1,j); xj2=lines(2,j);
        yj1=lines(3,j); yj2=lines(4,j);

        slopei=(yi2-yi1)/(xi2-xi1);
        slopej=(yj2-yj1)/(xj2-xj1);
        a=yi2-yi1;
        b=xi1-xi2;
        c=(xi2*yi1)-(xi1*yi2);
        dist1=abs(a*xj1+b*yj1+c)/sqrt(a^2+b^2);
        dist2=abs(a*xj2+b*yj2+c)/sqrt(a^2+b^2);
        if(isApproximatelyEqual(slopei,0,0.1)&&isApproximatelyEqual(slopej,0,0.1)&&isApproximatelyEqual(yi1,yj1,1)&&isApproximatelyEqual(yi2,yj2,1))
            if(xi1>xj1)
                %minx=xj1;
                newline=[xj1;xi2;yj1;yi2];
                mergelines=[mergelines,newline];
            else
                %minx=xi1;
                newline=[xi1;xj2;yi1;yj2];
                mergelines=[mergelines,newline];
            end
        elseif(isApproximatelyEqual(slopei,slopej,0.1)&&isApproximatelyEqual(dist1,0,1)&&isApproximatelyEqual(dist2,0,1))%另一条直线上两个点到第一条直线的距离近似0
            if(xi1>xj1)
                %minx=xj1;
                newline=[xj1;xi2;yj1;yi2];
                mergelines=[mergelines,newline];
            else
                %minx=xi1;
                newline=[xi1;xj2;yi1;yj2];
                mergelines=[mergelines,newline];
            end
            
        elseif(slopei==inf&&slopej==inf&&isApproximatelyEqual(xi1,xj1,1)&&isApproximatelyEqual(xi2,xj2,1))
                if(yi1>yj1)
                    %miny=yj1;
                    newline=[xj1;xi2;yj1;yi2];
                    mergelines=[mergelines,newline];
                else
                    %miny=yi1;
                    newline=[xi1;xj2;yi1;yj2];
                    mergelines=[mergelines,newline];
                end
            else
                newline=[xi1;xi2;yi1;yi2];
                mergelines=[mergelines,newline];
            end
        end
    end  
end




%         cosline1=[xi2-xi1,yj2-yi1];
%         cosline2=[xj1-xi1,yj1-yi1];
%         cosline3=[xj2-xi1,yj2-yi1];
%         cos1=dot(cosline1,cosline2)/(norm(cosline1)*norm(cosline2));
%         cos2=dot(cosline1,cosline3)/(norm(cosline1)*norm(cosline3));
%         if(isApproximatelyEqual(abs(cos1),1,0.01)&&isApproximatelyEqual(abs(cos2),1,0.01))
%             if(xi1<xj1)
%                 newline=[xi1;xj2;yi1;yj2];
%                 mergelines=[mergelines,newline];
%             elseif(xi1>xj1)
%                     newline=[xj1;xi2;yj1;yi2];
%                     mergelines=[mergelines,newline];
%                 else
%                     continue;
%             end
%         else
%             newline=[xi1;xi2;yj1;yj2];
%             mergelines=[mergelines,newline];
%         end
%     end
% end
% end

