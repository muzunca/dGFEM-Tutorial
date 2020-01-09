
%function [Diff_global,Conv_global,Reac_global,Fglobal]=global_system(mesh,fdiff,fadv,freact,fsource,
%                                                                     DBCexact,NBCexact,penalty,eps,degree)
%
% Purpose : Assembling local contributions to global matrices and right-hand
%          side
%
% Input :  
%         mesh     : mesh structure
%         fdiff    : diffusion function
%         fadv     : advection function
%         freact   : reaction function
%         fsource  : source function
%         DBCexact : Drichlet Boundary Condition
%         NBCexact : Neumann Boundary Condition
%         penalty  : penalty parameter to stabilize jumps over the edges
%            eps   : type of primal eps
%                    eps=1 ------>  NIPG
%                    eps=-1 ------> SIPG
%                    eps=0 -------> IIPG
%         degree   : degree of polynomials
%
% Output :
%         Diff_global : Nloc*Nel x Nloc*Nel
%                       global matrix for diffusion part
%         Conv_global : Nloc*Nel x Nloc*Nel
%                       global matrix for convection part            
%         Reac_global : Nloc*Nel x Nloc*Nel
%                       global matrix for reaction part
%         Fglobal     : Nloc*Nel x 1
%                       global vector for right-hand side


function [Diff_global,Conv_global,Reac_global,Fglobal]=global_system(mesh,fdiff,fadv,freact,fsource,...
                                                                DBCexact,NBCexact,penalty,eps,degree)

Nloc=(degree+1)*(degree+2)/2;            % local dimension
Nel=size(mesh.Elements,1);               % number of elements

%Initialize global matrix and right-hand side 
Diff_global=sparse(Nloc*Nel,Nloc*Nel);
Conv_global=sparse(Nloc*Nel,Nloc*Nel);
Reac_global=sparse(Nloc*Nel,Nloc*Nel);
Fglobal=sparse(Nloc*Nel,1);

% Assemble volume contributions

% Compute local volume matrix
[Diff_loc,Conv_loc,Reac_loc,Floc]=localmat_vol(mesh,fdiff,fadv,freact,fsource,degree);
  
Diff_loc=permute(Diff_loc,[3,2,1]);
Conv_loc=permute(Conv_loc,[3,2,1]);
Reac_loc=permute(Reac_loc,[3,2,1]);
Floc=permute(Floc,[3,2,1]);
   
kaux=0:Nloc:Nloc*(Nel-1);

id=reshape(repmat(1:Nloc,size(kaux',1),1),Nloc*size(kaux',1),1)+repmat(kaux',Nloc,1);
jd=ones(size(kaux',1)*Nloc,1);
F=Floc(:);

ie=[]; je=[]; D=[]; C=[]; 

% Interior Edges

% Extract  interior edges 
iedge=(mesh.intEdges);

%Get neighbors of edge
edge=mesh.EdgeEls(iedge,:);
E1=edge(:,1);
E2=edge(:,2);

kaux1=(E1-1)*Nloc;
kaux2=(E2-1)*Nloc;

% Compute local  matrices caused by interior edges
[B11,B22,B12,B21,C11,C22,C12,C21]= localmat_face(mesh,fdiff,fadv,penalty,eps,degree);

B11=permute(B11,[3,2,1]);
B22=permute(B22,[3,2,1]);
B12=permute(B12,[3,2,1]);
B21=permute(B21,[3,2,1]);

C11=permute(C11,[3,2,1]);
C22=permute(C22,[3,2,1]);
C12=permute(C12,[3,2,1]);
C21=permute(C21,[3,2,1]);

% Boundary edges :

% The penalty parameter is equal to 2*penal for both IIPG and SIPG and
% equal to penal for NIPG.

 switch eps
     case -1
       penalty=2*penalty;
     case 0
       penalty=2*penalty;
 end

if(~isempty(mesh.DbdEdges))
   % Extract Drichlet boundary edges 
   Dbedge=mesh.DbdEdges;

   %Get neighbors of edge
   edge=mesh.EdgeEls(Dbedge,:);
   E1=edge(:,1);

   kauxDB=(E1-1)*Nloc;

   % Compute local  matrices and right-hand side vector  caused by boundary edges
   [B11B,FlocDB,C11B]=localmat_Dbdyface(mesh,fdiff,fadv,DBCexact,penalty,eps,degree);

   B11B=permute(B11B,[3,2,1]);
   C11B=permute(C11B,[3,2,1]);
   FlocDB=permute(FlocDB,[3,2,1]);

   ie=[ie;reshape(repmat(1:Nloc,size(kauxDB,1)*Nloc,1),Nloc*Nloc*size(kauxDB,1),1)+repmat(kauxDB,Nloc*Nloc,1)];
   je=[je;repmat(reshape(repmat(1:Nloc,size(kauxDB,1),1),Nloc*size(kauxDB,1),1),Nloc,1)+repmat(kauxDB,Nloc*Nloc,1)];

   D=[D;B11B(:)];
   C=[C;C11B(:)];

   id=[id;reshape(repmat(1:Nloc,size(kauxDB,1),1),Nloc*size(kauxDB,1),1)+repmat(kauxDB,Nloc,1)];
   jd=[jd;ones(size(kauxDB,1)*Nloc,1)];

   F=[F;FlocDB(:)];
end

if(~isempty(mesh.NbdEdges))
   % Extract Neumann boundary edges 
   Nbedge=mesh.NbdEdges;

   %Get neighbors of edge
   edge=mesh.EdgeEls(Nbedge,:);
   E1=edge(:,1);

   kauxNB=(E1-1)*Nloc;

   % Compute local  matrices and right-hand side vector  caused by boundary edges
   FlocNB=localmat_Nbdyface(mesh,fdiff,NBCexact,degree);

   FlocNB=permute(FlocNB,[3,2,1]);

   id=[id;reshape(repmat(1:Nloc,length(kauxNB),1),Nloc*length(kauxNB),1)+repmat(kauxNB,Nloc,1)];

   jd=[jd;ones(length(kauxNB)*Nloc,1)];

   F=[F;FlocNB(:)];

end

ie=[reshape(repmat(1:Nloc,size(kaux',1)*Nloc,1),Nloc*Nloc*size(kaux',1),1)+repmat(kaux',Nloc*Nloc,1);
    reshape(repmat(1:Nloc,size(kaux1,1)*Nloc,1),Nloc*Nloc*size(kaux1,1),1)+repmat(kaux1,Nloc*Nloc,1);
    reshape(repmat(1:Nloc,size(kaux1,1)*Nloc,1),Nloc*Nloc*size(kaux1,1),1)+repmat(kaux2,Nloc*Nloc,1);
    reshape(repmat(1:Nloc,size(kaux1,1)*Nloc,1),Nloc*Nloc*size(kaux1,1),1)+repmat(kaux1,Nloc*Nloc,1);
    reshape(repmat(1:Nloc,size(kaux1,1)*Nloc,1),Nloc*Nloc*size(kaux1,1),1)+repmat(kaux2,Nloc*Nloc,1);ie];

 je=[repmat(reshape(repmat(1:Nloc,size(kaux',1),1),Nloc*size(kaux',1),1),Nloc,1)+repmat(kaux',Nloc*Nloc,1);
     repmat(reshape(repmat(1:Nloc,size(kaux1,1),1),Nloc*size(kaux1,1),1),Nloc,1)+repmat(kaux1,Nloc*Nloc,1);
     repmat(reshape(repmat(1:Nloc,size(kaux1,1),1),Nloc*size(kaux1,1),1),Nloc,1)+repmat(kaux2,Nloc*Nloc,1); 
     repmat(reshape(repmat(1:Nloc,size(kaux1,1),1),Nloc*size(kaux1,1),1),Nloc,1)+repmat(kaux2,Nloc*Nloc,1);  
     repmat(reshape(repmat(1:Nloc,size(kaux1,1),1),Nloc*size(kaux1,1),1),Nloc,1)+repmat(kaux1,Nloc*Nloc,1);je];
  
 ir=reshape(repmat(1:Nloc,size(kaux',1)*Nloc,1),Nloc*Nloc*size(kaux',1),1)+repmat(kaux',Nloc*Nloc,1);
 jr=repmat(reshape(repmat(1:Nloc,size(kaux',1),1),Nloc*size(kaux',1),1),Nloc,1)+repmat(kaux',Nloc*Nloc,1);
 
 
 D=[Diff_loc(:);B11(:);B22(:);B12(:);B21(:);D];
 C=[Conv_loc(:);C11(:);C22(:);C12(:);C21(:);C];
 
 Diff_global=Diff_global+sparse(ie,je,D,Nloc*Nel,Nloc*Nel);
 Conv_global=Conv_global+sparse(ie,je,C,Nloc*Nel,Nloc*Nel);
 Reac_global=Reac_global+sparse(ir,jr,Reac_loc(:),Nloc*Nel,Nloc*Nel);

 Fglobal=Fglobal+sparse(id,jd,F,Nloc*Nel,1);

return

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Dloc,Cloc,Rloc,Floc]=localmat_vol(mesh,fdiff,fadv,freact,fsource,degree)

% function [Dloc,Floc,Cloc,Rloc]=localmat_vol(mesh,fdiff,fadv,freact,fsource,degree)
%
% Purpose : Computes the local matrices caused by diffusion,convection and
%           reaction terms and the local right-hand side for each mesh element E
%
% Input   :
%         mesh     : mesh structure
%         fdiff    : diffusion function
%         fadv     : advection function
%         freact   : reaction function
%         fsource  : source function
%         degree   : degree of polynomials
%
% Output  :
%           Dloc : Nloc x Nloc x Nel
%                  local matrix for Diffision term
%           Cloc : Nloc x Nloc x Nel
%                  local matrix for convection term
%           Rloc : Nloc x Nloc x Nel
%                  local matrix for reaction term
%           Floc : Nloc x 1 x Nel
%                  local vector for right-hand side 

%Get quadrature points and weights on reference triangle
[nodes_ref,wgt]=quadrature(10);

Nloc=(degree+1)*(degree+2)/2;          % local dimension 
Nel=size(mesh.Elements,1);             % number of elements

%Compute values and derivatives of basis functions and determinant 
%and compute  global coordinates of quadarature point
[val_basis,der_basisx,der_basisy,determ,xx]=elem_basis(mesh,degree,nodes_ref);

%Quadrature weights
wgt=repmat(wgt,[1 1 Nel]);    % add Nel copies of the weights

% evaluate the diffusion function at the quadrature nodes
diff=feval(fdiff, xx(:,1,:),xx(:,2,:));
% evaluate the advection function at the quadrature nodes
[adv1,adv2]=feval(fadv, xx(:,1,:),xx(:,2,:));
% evaluate the reaction function at the quadrature nodes
reac =feval(freact,xx(:,1,:),xx(:,2,:));

% Compute  source function 
source=feval(fsource,fdiff,fadv,freact,xx(:,1,:),xx(:,2,:));

% Initialize to zero  the local matrix and local right-hand side 
Dloc=zeros(Nloc,Nloc,Nel);
Cloc=zeros(Nloc,Nloc,Nel);
Rloc=zeros(Nloc,Nloc,Nel);
Floc=zeros(Nloc,1,Nel);

% weights * determ and compute the transpose       
vol=permute(wgt.*determ,[2,1,3]);

clear wgt
   
%Compute  local matrices
for i=1:Nloc
    for j=1:Nloc
            
       % Diffusion part
       Dloc(i,j,:) =multiprod(vol,(der_basisx(:,j,:).*der_basisx(:,i,:)...
                                   +der_basisy(:,j,:).*der_basisy(:,i,:)).*diff);
                                
       % Convection part  
       Cloc(i,j,:) = multiprod(vol,(adv1.*der_basisx(:,j,:).*val_basis(:,i,:)...
                                   +adv2.*der_basisy(:,j,:).*val_basis(:,i,:)));
                                       
       % Reaction part 
       Rloc(i,j,:) =multiprod(vol,(val_basis(:,j,:).*val_basis(:,i,:)).*reac);
            
    end
        
    %Right-side
    Floc(i,1,:) =multiprod(vol,(val_basis(:,i,:).*source));
end

return;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D11,D22,D12,D21,C11,C22,C12,C21]=localmat_face(mesh,fdiff,fadv,penalty,eps,degree)

% function [D11,D22,D12,D21,C11,C22,C12,C21]=localmat_face(mesh,fdiff,fadv,penalty,eps,degree)
%
% Purpose : Compute local matrices and right-hand side  obtained by 
%           integration over interior  edges
%
% Input :  
%         mesh     : mesh structure
%         fdiff    : diffusion function
%         fadv     : advection function
%         penalty  : penalty parameter to stabilize jumps over the edges
%            eps   : type of primal eps
%                    eps=1 ------>  NIPG
%                    eps=-1 ------> SIPG
%                    eps=0 -------> IIPG
%          degree  : polynomila degree
%
% Output :
%         D11,D22,D12,D21  : Nloc x Nloc x Ned
%                            local matrix caused by diffusion part  
%         C11,C22,C12,C21  : Nloc x Nloc x Ned
%                            local matrix caused by convection part
%

global Equation

% Getr interior edges
iedge=(mesh.intEdges);

Nloc=(degree+1)*(degree+2)/2;         % local dimension
Ned=size(iedge,1);                    % number of interior edges

%Initialize the quadrature weights and points for edges
Nqu=12;  [nodes_ref,wge]=get_quadrature_segment(Nqu);

%Get neighbors of interior edges
edge=mesh.EdgeEls(iedge,:);
E1=edge(:,1);
E2=edge(:,2);

%Compute normal vector to edges E1
[normal_vec,area]=getNormal(mesh,E1,iedge);
normal_vec=permute(repmat(normal_vec,[1,1,Nqu]),[3,2,1]);
normal1=normal_vec(:,1,:);
normal2=normal_vec(:,2,:);

%Compute local coordinates of quadrature point on E1 and E2
s1=loc_coor_quad(mesh,iedge,E1,nodes_ref); 
s2=loc_coor_quad(mesh,iedge,E2,nodes_ref);

%Compute values and derivatives of basis functions and determinant 
%and compute  global coordinates of quadarature point on E1 on E2
[val_basis1,der_basis1x,der_basis1y,~,xx]=elem_basisf(mesh,degree,E1,s1);
[val_basis2,der_basis2x,der_basis2y,~,~]=elem_basisf(mesh,degree,E2,s2);

% Initialize to zero  the local matrices caused by diffusion
D11=zeros(Nloc,Nloc,Ned);
D22=zeros(Nloc,Nloc,Ned);
D12=zeros(Nloc,Nloc,Ned);
D21=zeros(Nloc,Nloc,Ned);

% Initialize to zero  the local matrices caused by convection part
C11=zeros(Nloc,Nloc,Ned);
C22=zeros(Nloc,Nloc,Ned);
C12=zeros(Nloc,Nloc,Ned);
C21=zeros(Nloc,Nloc,Ned);
       
% evaluate the diffusion function at the quadrature nodes
diff=feval(fdiff,xx(:,1,:),xx(:,2,:));
% evaluate the advection function at the quadrature nodes
[adv1,adv2]=feval(fadv,xx(:,1,:), xx(:,2,:));

%Penalty parameter
penalty=(penalty*ones(Ned,1))./((area.^Equation.b0).*ones(Ned,1));
penalty=permute(repmat(penalty,[1,1,Nqu]),[3,2,1]);
penalty=penalty.*diff;

%Quadrature weights
wge=repmat(wge,[1 1 Ned]);       %add Ned copies of the weights

area=permute(repmat(area,[1,1,Nqu]),[3,2,1]);
% weights * area and compute the transpose
area=permute(0.5*area.*wge,[2,1,3]);

%Define inflow and outflow edges
c=(adv1.*normal1+adv2.*normal2);
a=find(c(1,:,:)<0);
b=find(c(1,:,:)>=0);

 for i=1:Nloc
     for j=1:Nloc
      
     %Compute the entries of local matrix D11 
     
     T11=(normal1.*der_basis1x(:,j,:)+normal2.*der_basis1y(:,j,:)).*val_basis1(:,i,:).*(-0.5*diff)...
         +(normal1.*der_basis1x(:,i,:)+normal2.*der_basis1y(:,i,:)).*val_basis1(:,j,:).*(eps*(0.5)*diff)...
         +(penalty.*val_basis1(:,i,:).*val_basis1(:,j,:));
     
     D11(i,j,:) = multiprod(area,T11);
           
     %Compute the entries of local matrix D22 

     T22=(normal1.*der_basis2x(:,j,:)+normal2.*der_basis2y(:,j,:)).*val_basis2(:,i,:).*(0.5*diff)...
         +(normal1.*der_basis2x(:,i,:)+normal2.*der_basis2y(:,i,:)).*val_basis2(:,j,:).*(-0.5*eps*diff)...
         +(penalty.*val_basis2(:,i,:).*val_basis2(:,j,:));
     
     D22(i,j,:) = multiprod(area,T22);
          
         
     %Compute the entries of local matrix D12  

     T12=(normal1.*der_basis2x(:,j,:)+normal2.*der_basis2y(:,j,:)).*val_basis1(:,i,:).*(-0.5*diff)...
         +(normal1.*der_basis1x(:,i,:)+normal2.*der_basis1y(:,i,:)).*val_basis2(:,j,:).*(-0.5*eps*diff)...
         -(penalty.*val_basis1(:,i,:).*val_basis2(:,j,:));
     
     D12(i,j,:) = multiprod(area,T12);
   
     %Compute the entries of local matrix D21
     
     T21=(normal1.*der_basis1x(:,j,:)+normal2.*der_basis1y(:,j,:)).*val_basis2(:,i,:).*(0.5*diff)...
         +(normal1.*der_basis2x(:,i,:)+normal2.*der_basis2y(:,i,:)).*val_basis1(:,j,:).*(0.5*eps*diff)...
         -(penalty.*val_basis2(:,i,:).*val_basis1(:,j,:));
     
     D21(i,j,:) = multiprod(area,T21);

     
     % Compute convection part C11, C12, C22 and C21
     
     T1=(abs(c(:,:,a)).*val_basis1(:,j,a).*val_basis1(:,i,a));
     C11(i,j,a)=multiprod(area(:,:,a),T1);
        
     T2=(-abs(c(:,:,a)).*val_basis2(:,j,a).*val_basis1(:,i,a));
     C12(i,j,a)=multiprod(area(:,:,a),T2);
     
     T3=(-abs(c(:,:,b)).*val_basis1(:,j,b).*val_basis2(:,i,b));
     C21(i,j,b)=multiprod(area(:,:,b),T3);
     
     T4=(abs(c(:,:,b)).*val_basis2(:,j,b).*val_basis2(:,i,b));
     C22(i,j,b)=multiprod(area(:,:,b),T4);
     
     end
 end
              
return;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [D11,Floc,C11]=localmat_Dbdyface(mesh,fdiff,fadv,DBCexact,penalty,eps,degree)

% function [D11,Floc,C11]=localmat_Dbdyface(mesh,fdiff,fadv,DBCexact,penalty,eps,degree)
%
% Purpose : Compute local matrices and right-hand side  obtained by 
%           integration over Drichlet boundary edges
%
% Input :  
%         mesh     : mesh structure
%         fdiff    : diffusion function
%         fadv     : advection function
%         DBCexact : Drichlet Boundary Condition    
%         penalty  : penalty parameter to stabilize jumps over the edges
%            eps   : type of primal eps
%                    eps=1 ------>  NIPG
%                    eps=-1 ------> SIPG
%                    eps=0 -------> IIPG
%         degree   : degree of polynomials
%
% Output :
%         D11  : Nloc x Nloc x Ned
%                local matrix caused by diffusion part  
%         Floc : Nloc x 1 x Ned
%                local right-hand side vector
%         C11  : Nloc x Nloc x Ned
%                local matrix caused by convection part
%
global Equation;

Nloc=(degree+1)*(degree+2)/2;       % local dimension

% Initialize the quadrature weights and points for edges
Nqu=12;  [nodes_ref,wge]=get_quadrature_segment(Nqu);

Dbedge=mesh.DbdEdges;            % get Drichlet BC on edges
DNed=length(Dbedge);             % size of Drichlet BC

% Get neighbors of edge
edge=mesh.EdgeEls(Dbedge,:);
E1=edge(:,1);

% Compute normal vector to edges E1
[normal_vec,area]=getNormal(mesh,E1,Dbedge);
normal_vec=permute(repmat(normal_vec,[1,1,Nqu]),[3,2,1]);
normal1=normal_vec(:,1,:);
normal2=normal_vec(:,2,:);

% Compute local coordinates of quadrature points
s1=loc_coor_quad(mesh,Dbedge,E1,nodes_ref); 

% Compute values and derivatives of basis functions and determinant 
% and compute  global coordinates of quadarature point
[val_basis1,der_basis1x,der_basis1y,~,xx]=elem_basisf(mesh,degree,E1,s1);

% evaluate the diffusion function at the quadrature nodes
diff=feval(fdiff, xx(:,1,:),xx(:,2,:));
% evaluate the advection function at the quadrature nodes
[adv1,adv2]=feval(fadv, xx(:,1,:),xx(:,2,:));

% Penalty parameter
penalty=(penalty*ones(DNed,1))./((area.^Equation.b0).*ones(DNed,1));
penalty=permute(repmat(penalty,[1,1,Nqu]),[3,2,1]);
penalty=penalty.*diff;

% Quadrature weights
wge=repmat(wge,[1 1 DNed]);             % add DNed copies of the weights
area=permute(repmat(area,[1,1,Nqu]),[3,2,1]);
% weights * area and compute the transpose
area=permute(0.5*area.*wge,[2,1,3]);

% Initialize to zero  the local matrices
D11=zeros(Nloc,Nloc,DNed);
C11=zeros(Nloc,Nloc,DNed);
Floc=zeros(Nloc,1,DNed);
 
% Determine inflow edges
c=(adv1.*normal1+ adv2.*normal2);
a=find(c(1,:,:)<0);

% Compute  dirichlet BC
exactv=feval(DBCexact,fdiff,xx(:,1,:),xx(:,2,:));
 
for i=1:Nloc
    
   for j=1:Nloc
        
       %Compute the entries of local matrix D11   
       T11=(normal1.*der_basis1x(:,j,:).*val_basis1(:,i,:)+normal2.*der_basis1y(:,j,:).*val_basis1(:,i,:)).*(-diff)...
           +(normal1.*der_basis1x(:,i,:).*val_basis1(:,j,:)+normal2.*der_basis1y(:,i,:).*val_basis1(:,j,:)).*(eps*diff)...
           +(penalty.*val_basis1(:,i,:).*val_basis1(:,j,:));
        
       D11(i,j,:) = multiprod(area,T11);
                
       %Compute the entries of local matrix C11
       T=(abs(c(:,:,a)).*val_basis1(:,j,a).*val_basis1(:,i,a));
       C11(i,j,a) = multiprod(area(:,:,a),T);  
       
   end
     
    % Compute entries of local vector Floc  
     T1=(normal1.*der_basis1x(:,i,:).*exactv+normal2.*der_basis1y(:,i,:).*exactv).*(eps.*diff)...
        +(penalty.*val_basis1(:,i,:).*exactv);
    
     Floc(i,1,:)=multiprod(area,T1);
     Floc(i,1,a)=Floc(i,1,a)+multiprod(area(:,:,a),(abs(c(:,:,a)).*val_basis1(:,i,a).*exactv(:,:,a))); 
     
     
end
              
return;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Floc]=localmat_Nbdyface(mesh,fdiff,NBCexact,degree)

% function Floc=localmat_Nbdyface(mesh,fdiff,NBCexact,degree)
%
% Purpose : Compute local matrices and right-hand side  obtained by 
%           integration over Neumann boundary edges
%
% Input :  
%         mesh       : mesh structure
%         fdiff      : diffusion parameter
%         NBCexact   : Neumann Boundary condition     
%         degree     : degree of polynomials
%
% Output :  
%         Floc : Nloc x 1 x NNed
%                local right-hand side vector

Nloc=(degree+1)*(degree+2)/2;       % local dimension

% Initialize the quadrature weights and points for edges
Nqu=12;  [nodes_ref,wge]=get_quadrature_segment(Nqu);

Nbedge=mesh.NbdEdges;          % get Neumann BC on edges
NNed=length(Nbedge);           % size of Neumann BC

% Get neighbors of edge
edge=mesh.EdgeEls(Nbedge,:);
E1=edge(:,1);

% Compute local coordinates of quadrature points
s1=loc_coor_quad(mesh,Nbedge,E1,nodes_ref); 

% Compute values and derivatives of basis functions and determinant 
% and compute  global coordinates of quadarature point
[val_basis1,~,~,~,xx]=elem_basisf(mesh,degree,E1,s1);

[~,area]=getNormal(mesh,E1,Nbedge);

% Quadrature weights
wge=repmat(wge,[1 1 NNed]);             % add NNed copies of the weights
area=permute(repmat(area,[1,1,Nqu]),[3,2,1]);
% weights * area and compute the transpose
area=permute(0.5*area.*wge,[2,1,3]);

% Initialize the local vector to zero  
Floc=zeros(Nloc,1,NNed);
 
% Compute Neumann BC
exactv=feval(NBCexact,mesh,fdiff,xx(:,1,:),xx(:,2,:));
 
for i=1:Nloc
  % Compute entries of local vector Floc  
  Floc(i,1,:)=multiprod(area,val_basis1(:,i,:).*exactv);
end
              
return;

end






