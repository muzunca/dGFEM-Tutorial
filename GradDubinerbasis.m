

% function [dmodedr, dmodeds] = GradDubinerbasis(x,y,id,jd)
%
% Purpose: Return the derivatives of the Dubiner  basis (id,jd)
%          on the 2D simplex at (x,y). 

% P_(ij)=JacobiP(a,0,0,i)*(1-b)^(i) * JacobiP(b,2*i+1,0,j)
% 
% See "Analysis and Application of an Orthogonal Nodal Basis on
% Triangles for Discontinuous Spectral Element Methods", Shaozhong Deng
% and Wei Cai



function [dmodedr, dmodeds] = GradDubinerbasis(x,y,id,jd)

d=size(x);
a=((2*x)./(1-y))-1;   b=2*y-1;

fa =reshape( JacobiP(a, 0, 0, id),d); 
dfa =reshape( GradJacobiP(a, 0, 0, id),d);

gb = reshape(JacobiP(b, 2*id+1,0, jd),d);
dgb =reshape(GradJacobiP(b, 2*id+1,0, jd),d);

% x-derivative
% d/dx = da/dx d/da + db/dx d/db
dmodedr=dfa.*gb;

if(id>0)
  dmodedr = 2^(id+1).*dmodedr.*((1-y).^(id-1));
end


% y-derivative
% d/dy = da/dy d/da + db/dy d/db

if(id==0)
dmodeds=dfa.*gb.*(2*x)+fa.*dgb*2;
end

if(id>0)
 dmodeds =dfa.*gb.*(2*((1-b).^(id-1)).*(a+1));
 dmodeds =dmodeds+fa.*dgb.*(2^(id+1)*(1-y).^(id));
 dmodeds =dmodeds-fa.*gb.*((2^id)*id*(1-y).^(id-1));
end


return;
