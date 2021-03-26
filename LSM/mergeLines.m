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

function [ L ] = mergeLines( D,tau_theta,xi_s)

L=D;
L1=cell(1,2);
L2=cell(1,2);

[lengths,angles]=get_line_properties(D);
num_lines=size(D,1);

while 1
    clear prf;
    %sort lines in L in descending order of length (Psuedo Code line #4)
    [L, lengths, angles]=sort_lines_lengths(L,lengths,angles);
    
    %(x1,y1), (x2,y2) are the end points for each line segment in L
    x1=L(:,1); x2=L(:,3);
    y1=L(:,2); y2=L(:,4);
    
    %indL1 contains the index for line L1 and L1 is the line segment
    for indL1=1:size(L,1)
        if indL1>size(L,1)
            break;          %all lines have been processed
        end
        
        %exhaustively merge current line and then remove it from the active set
        angL1=angles(indL1);%its the angle for line segment L1
        
        %angular proximal lines (Psuedo Code line #8)
        inds=find(abs(angles-angL1)<tau_theta);%find lines with direction 'close to' angL1
        inds=[inds ; find(abs(angles-angL1)>(pi-2*pi/60))]; %to handle the wrap-around problem. (0,pi/60) should be compared with (pi-pi/60,pi)
        inds=setdiff(inds,indL1);%inds contains indicies for angulars proximal lines w.r.t L1
        if ~isempty(inds)
            cx1=x1(indL1); cx2=x2(indL1);
            cy1=y1(indL1); cy2=y2(indL1);
            L1{1}=[cx1 cy1];
            L1{2}=[cx2 cy2];
            
            %spatial proximal lines (Psuedo Code line #9)
            %from angular ones (inds contains indicies for angular proximal lines)
            %adaptive spatial proximal threshold
            tau_s=xi_s*lengths(indL1);
            x1diff=abs(repmat(cx1,length(inds),2)-[x1(inds) x2(inds)]);
            x2diff=abs(repmat(cx2,length(inds),2)-[x1(inds) x2(inds)]);
            inds1=inds(any([x1diff x2diff]<tau_s,2));
            if isempty(inds1)
                continue
            end
            y1diff=abs(repmat(cy1,length(inds1),2)-[y1(inds1) y2(inds1)]);
            y2diff=abs(repmat(cy2,length(inds1),2)-[y1(inds1) y2(inds1)]);
            inds2=inds1(any([y1diff y2diff]<tau_s,2));%inds2 contains angulars + spatial proximal lines
            
            % bad_inds contains the indicies to remove from L, (Psuedo Code line #10)
            bad_inds=[];
            
            for jj=1:length(inds2) %(inds2: spatial + angular proximal group indicies)
                %indL2 contains index for L2 and L2 is the line segment
                indL2=inds2(jj);
                L2{1}=[x1(indL2) y1(indL2)];
                L2{2}=[x2(indL2) y2(indL2)];
                
                %if mergeable, then merges line segments L1 and L2
                %as length and angle is available for both the lines, we pass that too (not part of psuedo code in paper)
                [M]=mergeTwoLines(L1,L2,xi_s,tau_theta,lengths([indL1,indL2]),angles([indL1,indL2]));
                
                if length(M)>1 %if M is not empty (Psuedo Code line #13)
                    %replace line L1 by the merged line M                 
                    L(indL1,:)=[M(1:4) mean(L([indL1,indL2],5))];
                    lengths(indL1)=M(5);
                    angles(indL1)=M(6);
                    L1{1}=M(1:2);
                    L1{2}=M(3:4);
       
                    %store indicies of all L2's to be removed from the set L
                    bad_inds=[bad_inds indL2]; % bad_inds contains the indicies to remove from set of lines L
                end
            end
            %update the set of lines L, x1, x2, y1, y2
            good_inds=setdiff(1:size(L,1),bad_inds); % removing indicies of merged lines, all L2 are mergered in L1 (Psuedo Code line #19)
            L=L(good_inds,:); % L now only containes the merged line, respective lines being merged are removed from L
            lengths=lengths(good_inds,:);
            angles=angles(good_inds,:);
            x1=L(:,1); x2=L(:,3);
            y1=L(:,2); y2=L(:,4);
        end
    end
    if size(L,1)==num_lines
        break;  %line merging has converged
    else
        num_lines=size(L,1);
    end
end
