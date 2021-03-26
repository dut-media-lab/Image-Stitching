function [inds]=getpoints(linenum,pointlist,except)

if~(exist('except','var'))
    except=0;
end

inds=[];
len=length(pointlist);
for i=1:len
    if( (~isempty(find(pointlist(i).lines==linenum, 1)))&& ...
            (isempty(find(pointlist(i).lines==except, 1))))
        inds=[inds i];
        
    end
end

end