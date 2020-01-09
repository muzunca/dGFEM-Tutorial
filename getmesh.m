
% function mesh = getmesh(Nodes,Elements,Dirichlet,Neumann,periodic)
%
% Purpose       : Generates initial mesh data structure.
% 
% INPUT 
%       Nodes     :  array of Nodes 
%       Elements  :  array of Element indices 
%       Dirichlet :  Dirichlet boundary edges 
%       Neumann   :  Neumann boundary edges 
%       periodic  :  Periodicity of boundary
%
% OUTPUT
%       mesh      :  current mesh
%
%  This routine is a modified firm of AFEM by L. Chen & C. Zhang 


function mesh = getmesh(Nodes,Elements,Dirichlet,Neumann)


% Label the mesh by the longest edge rule
Elements = label(Nodes,Elements);

% Define Edges, Interior Edges and Boundary Edges
totalEdge = sort([Elements(:,[2,3]); Elements(:,[3,1]); Elements(:,[1,2])],2);
[i,j,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
Edges = [j,i];            % Edges
bdyEdges = find(s==1);    % Boundary Edges
intEdges = find(s==2);      % Interior Edges

DbdEdges=[];
NbdEdges=[];

for j=1:size(Dirichlet,1)
 DbdEdges(j)=find(ismember(Edges,Dirichlet(j,:),'rows'));
end
for j=1:size(Neumann,1)
 NbdEdges(j)=find(ismember(Edges,Neumann(j,:),'rows'));
end   


% Vertices
vertices1=[Nodes(Elements(:,1),1),Nodes(Elements(:,1),2)];
vertices2=[Nodes(Elements(:,2),1),Nodes(Elements(:,2),2)];
vertices3=[Nodes(Elements(:,3),1),Nodes(Elements(:,3),2)];

vers=version('-release');
vers=uint8(vers(4));

% Define the elements containing the edge 
NT=size(Elements,1);
if vers>49
    [~, i2, j] = unique(totalEdge,'rows','legacy');
else
    [~, i2, j] = unique(totalEdge,'rows');
end
i1(j(3*NT:-1:1)) = 3*NT:-1:1; i1=i1';
k1 = ceil(i1/NT); t1 = i1 - NT*(k1-1);
k2 = ceil(i2/NT); t2 = i2 - NT*(k2-1);
EdgeEls = [t1,t2];

% Define the triangle using edges 
ElementsE = reshape(j,NT,3);  

% Generate mesh data structure
mesh = struct('Nodes',Nodes, 'Elements',Elements, 'Edges',Edges, 'DbdEdges',DbdEdges,'NbdEdges',NbdEdges,...
              'intEdges',intEdges,'vertices1', vertices1, 'vertices2' ,vertices2, 'vertices3' ,vertices3, ...
              'Dirichlet',Dirichlet, 'Neumann',Neumann,'EdgeEls',EdgeEls,'ElementsE',ElementsE,'bdyEdges',bdyEdges  );           
return;
      


