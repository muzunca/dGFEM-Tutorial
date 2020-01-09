
% function  s=loc_coor_quad(mesh,face,E,xg)
%
% Purpose : Compute local coordinates of the given quadrature points on
%           [-1,1] on edge of reference triangle
%
% Input  :
%          mesh : mesh structure
%          face : list of edges
%          E    : list of triangles 
%          xg   : quadrature points on [-1,1] segment
%
%
% Output :
%          s   : Nqu x 2 x Ned
%                local coordinates of the  quadrature points


function  s=loc_coor_quad(mesh,face,E,xg)

Ned=size(E,1);          % number of triangles in E
Nqu=size(xg,1);         % number of quadrature points

% Compute coordinates of each edge 
ax=mesh.Nodes(mesh.Edges(face,1),1);
ay=mesh.Nodes(mesh.Edges(face,1),2);

bx=mesh.Nodes(mesh.Edges(face,2),1);
by=mesh.Nodes(mesh.Edges(face,2),2);

x1=repmat(reshape((bx-ax)',1,1,Ned),[Nqu,1,1]);
x2=repmat(reshape((bx+ax)',1,1,Ned),[Nqu,1,1]);

y1=repmat(reshape((by-ay)',1,1,Ned),[Nqu,1,1]);
y2=repmat(reshape((by+ay)',1,1,Ned),[Nqu,1,1]);

xg=repmat(reshape(xg',1,1,Nqu),[Ned,1,1]);
xg=permute(xg,[3,2,1]);

% Firstly, obtain the values of quadrature points on physical edges 
s1=(xg.*x1+x2)*0.5;     % x-axis coordinate 
s2=(xg.*y1+y2)*0.5;     % y-axis coordinate 


% The transformation matrix
BE=[mesh.vertices2(E,:)-mesh.vertices1(E,:)  mesh.vertices3(E,:)-mesh.vertices1(E,:)];

% Compute determinat of transformation matrix BE
determ=BE(:,1).*BE(:,4)-BE(:,3).*BE(:,2);

% Binv
Binv  = [BE(:,4) -BE(:,2) -BE(:,3) BE(:,1)];
Binv  = Binv./repmat(determ,[1 size(Binv,2)]);

% Binv^T
BinvT = [Binv(:,1) Binv(:,3) Binv(:,2) Binv(:,4)];

% Now, compute local coordinates of quadrature points by using
% transformation matrix

bb1=s1-reshape(repmat(mesh.vertices1(E,1)',[Nqu 1]),Nqu,1,[]);
bb2=s2-reshape(repmat(mesh.vertices1(E,2)',[Nqu 1]),Nqu,1,[]);

% local x-axis coordinate 
ss11= bb1.*repmat(reshape(BinvT(:,1)',1,1,Ned),[Nqu,1,1])...
         +bb2.*repmat(reshape(BinvT(:,2)',1,1,Ned),[Nqu,1,1]);
 
% local y-axis coordinate 
ss12= bb1.*repmat(reshape(BinvT(:,3)',1,1,Ned),[Nqu,1,1])...
         +bb2.*repmat(reshape(BinvT(:,4)',1,1,Ned),[Nqu,1,1]);      
     
s=[ss11, ss12];

return ;

