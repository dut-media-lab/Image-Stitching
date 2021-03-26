function cp=crossproduct(p1,p2,o)
% if~(exist('o','var'))
%     o=[0 0];
% end
cp=(p1(1)-o(1))*(p2(2)-o(2))-(p2(1)-o(1))*(p1(2)-o(2));

end