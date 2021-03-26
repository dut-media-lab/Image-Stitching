function [longline]=houghline(imgpath)
I=imread(imgpath);    
I=rgb2gray(I);    
thresh=[0.01,0.17];   
sigma=2;   
f = edge(double(I),'canny',thresh,sigma);    
% figure(1),imshow(f,[]);    
% title('canny Detect Result');    
    
[H, theta, rho]= hough(f,'RhoResolution', 0.5);    
peak=houghpeaks(H,5);    
% hold on      
lines=houghlines(f,theta,rho,peak); 

figure,imshow(f,[]),title('Hough Transform Detect Result'),hold on    
for k=1:length(lines)    
    xy=[lines(k).point1;lines(k).point2];    
    plot(xy(:,1),xy(:,2),'LineWidth',4,'Color',[0.6, 0.6, 0.6]);    
end 

point1=cat(1,lines.point1);
point2=cat(1,lines.point2);
longline=[point1,point2];
longline=longline';
longline([2,3],:)=longline([3,2],:);



