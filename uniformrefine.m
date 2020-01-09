
% function [mesh] = uniformrefine(mesh)
%
% Purpose : refines the current triangulation by dividing
%           each triangle into four similar triangles
%
% Input:
%         mesh :  current mesh
% 
% Output: 
%         mesh :  new mesh after refinement
%
%  This routine is a modified form of AFEM by L. Chen & C. Zhang 

function [mesh] = uniformrefine(mesh)


% Construct data structure
edge = [mesh.Elements(:,[1,2]); mesh.Elements(:,[1,3]); mesh.Elements(:,[2,3])];
edge = unique(sort(edge,2),'rows');
N = size(mesh.Nodes,1); NT = size(mesh.Elements,1); NE = size(edge,1);
d2p = sparse(edge(:,[1,2]),edge(:,[2,1]),[1:NE,1:NE],N,N);


% New Nodess from the mid points of each edge
newNodes = (mesh.Nodes(edge(:,1),:)+mesh.Nodes(edge(:,2),:))/2; 
mesh.Nodes = [mesh.Nodes; newNodes]; 
marker = N+1:N+NE; 

%--------------------------------------------------------------------------
% refine each triangle into four triangles
%     3
%    / \
%   6 - 5
%  / \ / \
% 1 - 4 - 2
%--------------------------------------------------------------------------
t=1:NT;
p(t,1:3) = mesh.Elements(t,[1 2 3]); 
p(t,4) = marker(d2p(sub2ind(size(d2p),mesh.Elements(t,1), mesh.Elements(t,2))));
p(t,5) = marker(d2p(sub2ind(size(d2p),mesh.Elements(t,2), mesh.Elements(t,3))));
p(t,6) = marker(d2p(sub2ind(size(d2p),mesh.Elements(t,3), mesh.Elements(t,1))));

t=1:NT;
mesh.Elements(t,[1 2 3]) = [p(t,1) p(t,4) p(t,6)];
mesh.Elements(NT+1:2*NT,:) = [p(t,4) p(t,2) p(t,5)];
mesh.Elements(2*NT+1:3*NT,:) = [p(t,6) p(t,5) p(t,3)];
mesh.Elements(3*NT+1:4*NT,:) = [p(t,4) p(t,5) p(t,6)];



%--------------------------------------------------------------------------
% Update boundary edges
%--------------------------------------------------------------------------
mesh.Dirichlet = sort(updatebd(mesh.Dirichlet,marker,d2p),2);
mesh.Neumann = sort(updatebd(mesh.Neumann,marker,d2p),2);

mesh = getmesh(mesh.Nodes,mesh.Elements,mesh.Dirichlet,mesh.Neumann);

end

function bdEdge = updatebd(bdEdge,marker,d2p)
% UPDATEDBD refine the boundary edges
%
% INPUT
%   bdEdge:  set of boundary edges
%   marker:  new Nodes index for marked edge
%      d2p:  index mapping from dual edge to primary edge
% 
% OUTPUT
%   bdEdge:  set of refined boundary edges
%

NB = size(bdEdge,1);
for k = 1:NB
     i = bdEdge(k,1); 
     j = bdEdge(k,2);
     bdEdge(k,:) = [i,marker(d2p(i,j))];
     bdEdge(size(bdEdge,1)+1,:) = [marker(d2p(i,j)),j];
end

end



