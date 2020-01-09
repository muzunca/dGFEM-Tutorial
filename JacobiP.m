
%
% function [P] = JacobiP(x,alpha,beta,N)
%
% Purpose: Evaluate Jacobi Polynomial of type (alpha,beta) > -1
%          (alpha+beta <> -1) at points x for order N and 
%           returns P[1:length(xp))]
% Note   : They are normalized to be orthonormal.
%
%   This routine is part of the MATLAB  code that
%   accompanies "Nodal Discontinuous Galerkin Methods Algorithms, 
%   Analysis, and Applications" by Jan S. Hesthaven and Tim Warburton.



function [P] = JacobiP(x,alpha,beta,N)

% Turn points into row if needed.
xp = x(:) ;
dims = size(xp);
if (dims(2)==1)
    xp = xp'; 
end

PL = zeros(N+1,length(xp)); 

% Initial values P_0(x) and P_1(x)
PL(1,:) = 1;
if (N==0)
    P=PL';
   return;
end;
PL(2,:) = 0.5*(alpha-beta+(alpha+beta+2)*xp);
if (N==1) 
    P=PL(N+1,:)'; 
    return;
end

% Forward recurrence using the symmetry of the recurrence.
for i=2:N
h=2*(i-1)+alpha+beta;
h1=alpha+beta+(i-1)+1;
an=((h+1)*(h+2))/(2*(i-1+1)*h1);
bn=((alpha^2-beta^2)*(h+1))/(2*(i-1+1)*h*h1);
cn=((i-1+alpha)*(i-1+beta)*(h+2))/((i-1+1)*h1*h);

PL(i+1,:) = (an*xp+bn).*PL(i,:)-cn*PL(i-1,:);
end;

P =PL(N+1,:)';

return;