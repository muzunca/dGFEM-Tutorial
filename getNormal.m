
% [normal,area]=getNormal(mesh,Element,Edge)
%
% Purpose : compute the outward-pointing unit normal to Edge(i) edge, 
%           which belongs to the Ele(i) triangle of mesh 
%
% Input   :
%          mesh
%          Element   : list of triangles
%          Edge      : list of edges
% 
% Output  :
%         normal : the outward-pointing unit normal
%         area   : length of edge
%
%   This routine is part of the MATLAB Fem code that
%   accompanies "Understanding and Implementing the Finite
%   Element Method" by Mark S. Gockenbach (copyright SIAM 2006).

function [normal,area]=getNormal(mesh,Element,Edge)

j=zeros(size(Element,1),1);
e1=zeros(size(Element,1),1);
area=zeros(size(Element,1),1);
normal=zeros(size(Element,1),2);

for k=1:size(Element,1)
    j(k)=find(Edge(k)==abs(mesh.ElementsE(Element(k),:)));
    e1(k)=mesh.ElementsE(Element(k),j(k));
    
    if (j(k)+1>3)
        x=mod(j(k)+1,3);
    else
        x=j(k)+1;
    end
    if (j(k)+2>3)
        y=mod(j(k)+2,3);
    else
        y=j(k)+2;
    end
    if (mesh.Elements(Element(k),x)>mesh.Elements(Element(k),y))
    e1(k)=-e1(k);
    end
end

c=mesh.Nodes(mesh.Edges(Edge,1:2)',:);

a=find(e1>0);
b=find(e1<0);

if find(e1>0)
    n(a,:)=[c(2*(a(:)-1)+2,2)-c(2*(a(:)-1)+1,2), c(2*(a(:)-1)+1,1)-c(2*(a(:)-1)+2,1)];
end

if find(e1<0)
    n(b,:)=[c(2*(b(:)-1)+1,2)-c(2*(b(:)-1)+2,2), c(2*(b(:)-1)+2,1)-c(2*(b(:)-1)+1,1)];
end

for k=1:size(Element,1)
   area(k)=norm(n(k,:));
   normal(k,:)=n(k,:)/area(k);
end

return;

