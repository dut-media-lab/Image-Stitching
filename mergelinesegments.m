function [mergelines] = mergelinesegments(lines)
% merge line segments in image coordinates,
% lines: 6*M matrix M lines
% |x1,x2,y1,y2|;...

judge = zeros(1,size(lines,2));
lines = [lines;judge];
mergelines = [];

i = 1;
origin_num = size(lines,2);
num_lines = size(lines,2);

while i < num_lines
    for j=1:size(lines,2)
        if j == i || lines(7,j) == 1
            continue;
        else
            if lines(1,i) < lines(2,i)
                x1_i = lines(1,i);
                x2_i = lines(2,i);
                y1_i = lines(3,i);
                y2_i = lines(4,i);
                p1_i = [x1_i, y1_i];
                p2_i = [x2_i, y2_i];
            else
                x1_i = lines(2,i);
                x2_i = lines(1,i);
                y1_i = lines(4,i);
                y2_i = lines(3,i);
                p1_i = [x1_i, y1_i];
                p2_i = [x2_i, y2_i];
            end
            if lines(1,j) < lines(2,j)
                x1_j = lines(1,j);
                x2_j = lines(2,j);
                y1_j = lines(3,j);
                y2_j = lines(4,j);
                p1_j = [x1_j, y1_j];
                p2_j = [x2_j, y2_j];
            else
                x1_j = lines(2,j);
                x2_j = lines(1,j);
                y1_j = lines(4,j);
                y2_j = lines(3,j);
                p1_j = [x1_j, y1_j];
                p2_j = [x2_j, y2_j];
            end
            
            if x1_i==x2_i
                slope_i = Inf;
            else
                slope_i = (y2_i-y1_i)/(x2_i-x1_i);
            end
            if x1_j==x2_j
                slope_j = Inf;
            else
                slope_j = (y2_j-y1_j)/(x2_j-x1_j);
            end
            if slope_i ~= Inf && slope_j ~= Inf
                theta_ij = atand(abs((slope_j-slope_i)/(1+slope_i*slope_j)));
            elseif slope_i == Inf && slope_j ~= Inf
                theta_ij = 90-atand(slope_j);
            elseif slope_i ~= Inf && slope_j == Inf
                theta_ij = 90-atand(slope_i);
            else
                theta_ij = 0;
            end
            
            a_i = y2_i-y1_i;
            b_i = x1_i-x2_i;
            c_i = (x2_i*y1_i)-(x1_i*y2_i);
            dist1_i_j1 = abs(a_i*x1_j+b_i*y1_j+c_i)/sqrt(a_i^2+b_i^2);
            dist2_i_j2 = abs(a_i*x2_j+b_i*y2_j+c_i)/sqrt(a_i^2+b_i^2);
            distl_p1_i_p2_j = norm(p2_j-p1_i);
            dists_p1_j_p2_i = norm(p2_i-p1_j);
            if dists_p1_j_p2_i>distl_p1_i_p2_j
                temp = distl_p1_i_p2_j;
                distl_p1_i_p2_j = dists_p1_j_p2_i;
                dists_p1_j_p2_i = temp;
            end
            
            
            if( theta_ij<8 && dist1_i_j1<2 && dist2_i_j2<2 && distl_p1_i_p2_j<100*dists_p1_j_p2_i )
                if( x1_i>=x1_j && norm(p1_i-p2_j)<=0.5*norm(p2_i-p1_j) )
                    newline=[x1_j;x2_i;y1_j;y2_i;1;sqrt((x2_i-x1_j)^2+(y2_i-y1_j)^2);0];
                    mergelines = [mergelines,newline];
                    lines = [lines,newline];
                    lines(7,i) = 1;
                    lines(7,j) = 1;
                    if j > origin_num
                        mergelines(7,j-origin_num) = 1;
                    end
                elseif( x1_j>x1_i && norm(p1_j-p2_i)<=0.5*norm(p2_j-p1_i) )
                    newline=[x1_i;x2_j;y1_i;y2_j;1;sqrt((x2_j-x1_i)^2+(y2_j-y1_i)^2);0];
                    mergelines = [mergelines,newline];
                    lines = [lines,newline];
                    lines(7,i) = 1;
                    lines(7,j) = 1;
                    if j > origin_num
                        mergelines(7,j-origin_num) = 1;
                    end
                end
            end
        end
    end
    num_lines = size(lines,2);
    i = i+1;
end

if size(mergelines,1)>0
    notextra = mergelines(7,:)==0;
    mergelines = mergelines(1:6,notextra);
end
end

