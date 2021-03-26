function [newlines, T] = normaliselines(lines)

    %line translate to cordinate center
    delX=sum(lines(:,1))/sum(lines(:,3));
    delY=sum(lines(:,2))/sum(lines(:,3));
    trans = [1 0 -delX; 0 1 -delY; 0 0 1]; % 
    
    res=lines*trans';
    

  % scale transform  
  colsum=sum(res.*res,1);
  s=sqrt((colsum(1)+colsum(2))/2*colsum(3));
  
    
    T = [1       0       -delX
         0       1       -delY
         0       0         s   ];
    
    newlines = lines*T';
    
end