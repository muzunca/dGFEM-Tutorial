

% function [P] = Dubinerbasis(x,y,i,j)
% Purpose : Evaluate 2D orthogonal polynomial
%           on simplex at (x,y) of order (i,j).
%
% P_(ij)=JacobiP(a,0,0,i)*(1-b)^(i) * JacobiP(b,2*i+1,0,j)
% 
% See "Analysis and Application of an Orthogonal Nodal Basis on
% Triangles for Discontinuous Spectral Element Methods", Shaozhong Deng
% and Wei Cai


function [P] = Dubinerbasis(x,y,i,j)

d=size(x);
a=((2*x)./(1-y))-1;
b=2*y-1;

h1 = JacobiP(a,0,0,i);      h1=reshape(h1,d);
h2 = JacobiP(b,2*i+1,0,j);  h2=reshape(h2,d);
P = (2^i)*h1.*h2.*(1-y).^i;

return;
