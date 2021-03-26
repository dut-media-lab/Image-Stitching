function [] = meshdraw(warped_img, wX, wY, off, C1, C2)

figure;
imshow(warped_img);
hold on
if off(1)<1
    wX_draw = wX+off(1);
else
    off_x = min(wX(:,1));
    wX_draw = wX-off_x;
end
if off(2)<1
    wY_draw = wY+off(2);
else
    off_y = min(wY(1,:));
    wY_draw = wY-off_y;    
end
for i=1:C1+1
    for j=1:C2
        wX_draw_temp = [wX_draw(i,j); wX_draw(i,j+1)];
        wY_draw_temp = [wY_draw(i,j); wY_draw(i,j+1)];
        line(wX_draw_temp,wY_draw_temp,'LineWidth',1.5,'Color','Red');
        hold on;
    end
end
for i=1:C2+1
    for j=1:C1
        wX_draw_temp = [wX_draw(j,i); wX_draw(j+1,i)];
        wY_draw_temp = [wY_draw(j,i); wY_draw(j+1,i)];
        line(wX_draw_temp,wY_draw_temp,'LineWidth',1.5,'Color','Red');
        hold on;
    end
end
hold off

end