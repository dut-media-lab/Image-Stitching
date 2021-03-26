% This file is part of LSM.
% 
% LSM is free software: you can redistribute it and/or modify
% it under the terms of the GNU v3.0 General Public License as published by
% the Free Software Foundation.
% 
% LSM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details
% <http://www.gnu.org/licenses/>.

function [M]=mergeTwoLines(L1,L2,xi_s,tau_theta,lengths,angles)
%give more importance to the longer of the two lines
%longer line is named line L1 and angles are compared w.r.t L1
if lengths(2)>lengths(1)
    %swap lines
    tmp=L1;
    L1=L2;
    L2=tmp;
    %swap lengths
    tmp=lengths(2);
    lengths(2)=lengths(1);
    lengths(1)=tmp;
    %swap angels
    tmp=angles(2);
    angles(2)=angles(1);
    angles(1)=tmp;
end
p1L1=L1{1}; p2L1=L1{2}; % first and second point of line segment L1
p1L2=L2{1}; p2L2=L2{2}; % first and second point of line segment L2
l1=lengths(1); % length of line segment L1
l2=lengths(2); % length of line segment L2
theta1=angles(1); % angle of line segment L1
theta2=angles(2); % angle of line segment L2

%(Psuedo Code line #10)
dmat=[norm(p1L1-p1L2) norm(p1L1-p2L2) ; norm(p2L1-p1L2) norm(p2L1-p2L2)];   %pairwise distances between the 4 end-points
[d,i,j]=max2d(-dmat);
d=-d;
c1=L1{i(1)}; c2=L2{j(1)};           %closest end-points of lines L1 and L2 (c1,c2 as in algo)
L1o=L1{2/i(1)}; L2o=L2{2/j(1)};       %other end points of L1 and L2

tau_s=xi_s*l1; 
if d>tau_s
    M=-1;
else
    %computing TAU_THETA_STAR (Psuedo Code line #16)
    lambda=((l2/l1)+(d/tau_s));
    %2--->skewness of logistic sigmoid
    tau_theta_star= (1-(1./(1+exp(-2*(lambda-1.5)))))*tau_theta; %1.5--->1.5 out of 2 is assigned greater threshold.
    theta=abs(theta1-theta2);
    if theta < tau_theta_star || theta > pi-tau_theta_star
        prf_flag=0;
        %computing farthest end points (Psuedo Code line #19)
        dist = [pdist2(c1,L1o) pdist2(c1,c2) pdist2(c1,L2o) pdist2(L1o,L2o) pdist2(L1o,c2) pdist2(c2,L2o) ];
        [max_dist, max_ind]=max(dist);
        if max_ind==1
            M=[c1 L1o];
        elseif max_ind==2
            %can never happen
            M=[c1 c2];
        elseif max_ind==3
            prf_flag=1;
            M=[c1 L2o];
        elseif max_ind==4
            prf_flag=1;
            M=[L1o L2o];
        elseif max_ind==5
            prf_flag=1;
            M=[L1o c2];
        else
            %can never happen
            M=[c2 L2o];
        end
        
        if M~=-1           
            %gradient of merged line should be similar to gradient of line L1
            %final check on the absolute angular difference of the longer segment L1 and the merged segment M.
            %(Psuedo Code line #22-24)
            [gm m_dx m_dy]=line_gradient(M(1:2),M(3:4));
            thetaM=atan2(m_dy, m_dx);
            lengthM=sqrt(m_dx.^2+m_dy.^2);
            
            if thetaM<=0
                thetaM=thetaM+pi;
            end
             
            if abs(theta1-thetaM)<(tau_theta/2) 
                M(5)=lengthM;
                M(6)=thetaM;
            else
                M=-1;
            end
        end
    else
        M=-1;
    end
end

end

function [g, dx, dy]=line_gradient(p1,p2)
dx=p2(1)-p1(1);
if dx
    dy=p2(2)-p1(2);
    g=dy/dx;
else
    dy=p2(2)-p1(2);
    g=1e10;
end
end