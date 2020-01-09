%function [R,RJac]=nonlinear_global(coef,mesh,freact_nonlinear,degree)
%
% Purpose: Compute the vector for the integration to the nonlinear 
%          reaction and its Jacobian matrix at the current iterate
%          of coefficient vector "coef"
%
% Input:
%
%                coef  : current coefficient vector
%                mesh  : mesh structure
%    freact_nonlinear  : nonlinear reaction function
%              degree  : polynomial degree
%
% Output:
%
%       R   :   Nloc*Nel x 1
%               global nonlinear reaction vector
%
%    RJac   :   Nloc*Nel x Nloc*Nel
%               global Jacobian matrix to nonlinear reaction
%

function [R,RJac]=nonlinear_global(coef,mesh,freact_nonlinear,degree)

Nloc=(degree+1)*(degree+2)/2;   % local dimension
Nel=size(mesh.Elements,1);      % number of elements

%Initialize global R and RJac
RJac=sparse(Nloc*Nel,Nloc*Nel);
R=zeros(Nloc*Nel,1);

kaux=0:Nloc:Nloc*(Nel-1);
ir=reshape(repmat(1:Nloc,size(kaux',1)*Nloc,1),Nloc*Nloc*size(kaux',1),1)+repmat(kaux',Nloc*Nloc,1);
jr=repmat(reshape(repmat(1:Nloc,size(kaux',1),1),Nloc*size(kaux',1),1),Nloc,1)+repmat(kaux',Nloc*Nloc,1);
id=reshape(repmat(1:Nloc,size(kaux',1),1),Nloc*size(kaux',1),1)+repmat(kaux',Nloc,1);
jd=ones(size(kaux',1)*Nloc,1);

% Compute local volume matrices
[R_loc,RJac_loc]=nonlinear_vol(coef,mesh,freact_nonlinear,degree);
RJac_loc=permute(RJac_loc,[3 2 1]);
R_loc=permute(R_loc,[3 2 1]);

RJac=RJac+sparse(ir,jr,RJac_loc(:),Nloc*Nel,Nloc*Nel);
R=R+sparse(id,jd,R_loc(:),Nloc*Nel,1); 

return;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rloc,RJacloc]=nonlinear_vol(coef,mesh,freact_nonlinear,degree)

%Get quadrature points and weights on reference triangle
[nodes_ref,wgt]=quadrature(10);

Nloc=(degree+1)*(degree+2)/2;   % local dimension 
Nel=size(mesh.Elements,1);      % number of elements

%Compute values and derivatives of basis functions and determinant 
%and compute  global coordinates of quadarature point
[val_basis,~,~,determ,~]=elem_basis(mesh,degree,nodes_ref);

%Quadrature weights
wgt=repmat(wgt,[1 1 Nel]);    % add Nel copies of the weights

% Initialize to zero  the local matrix and local right-hand side 
RJacloc=zeros(Nloc,Nloc,Nel);
Rloc=zeros(Nloc,1,Nel);

% weights * determ and compute the transpose       
vol=permute(wgt.*determ,[2,1,3]);

% Reshape the coefficient vector
coef=reshape(coef,[Nloc 1 Nel]);

% Value of the solution for the current coefficient vector "coef"
val=multiprod(val_basis,coef);

% Compute the values of the nonlinear term and its derivative
[r,dr]=freact_nonlinear(val);

%Compute  local matrices  
 for i=1:Nloc
     for j=1:Nloc
         RJacloc(i,j,:) =multiprod(vol,val_basis(:,i,:).*val_basis(:,j,:).*dr);
     end

     Rloc(i,1,:) =multiprod(vol,val_basis(:,i,:).*r);
 end
  
    
  return;
end