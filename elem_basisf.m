
%function [vbasis,Pdbasisx,Pdbasisy,determ,xx]=elem_basisf(mesh,degree,E,nodes_ref)
%
% Purpose: Compute the values of the basis functions, 
%          global derivative of the basis functions,
%          the determinant of the transformation matrix and
%          the values of quadrature points on all edges of triangles.
%
% Input:
%         mesh      : mesh structure
%         degree    : polynomial degree
%         E         : Edges (it can be interior or boundary edges)
%         nodes_ref : Nqu x 2 x Ned
%                     quadrature points on all edges of reference
%                     triangles
%
% Output:
%         vbasis    : Nqu x Nloc x Ned
%                     values of Nloc basis function at the Nqu integration
%                     points for edges 1:Ned
%         Pdbasisx  : Nqu x Nloc x Ned
%                     x derivative of Nloc basis function at the Nqu
%                     quadrature points for edges 1:Ned
%         Pdbasisy  : Nqu x Nloc x Ned
%                     y derivative of Nloc basis function at the Nqu
%                     quadrature points for edges 1:Ned
%         determ    : Nqu x 1 x Ned
%                     determinant of the transformation matrix between
%                     reference triangle and physical triangle
%         xx        : Nqu x 2 x Ned
%                     the values of quadrature points on the physical
%                     edges



function [vbasis,Pdbasisx,Pdbasisy,determ,xx]=elem_basisf(mesh,degree,E,nodes_ref)

global Equation;

Ned=size(E,1);                    % number of edges
Nqu=size(nodes_ref,1);            % number of quadrature points
Nloc=(degree+1)*(degree+2)*0.5;   % local dimension


% Construct  the transformation matrix that maps from the reference triangle
% defined by (0,0),(0,1), (1,0) onto the triangle in the physical domain

% The transformation matrix
BE=[mesh.vertices2(E,:)-mesh.vertices1(E,:)  mesh.vertices3(E,:)-mesh.vertices1(E,:)];

% Compute determinat of transformation matrix BE
determ=BE(:,1).*BE(:,4)-BE(:,3).*BE(:,2);

% Binv
Binv  = [BE(:,4) -BE(:,2) -BE(:,3) BE(:,1)];
Binv  = Binv./repmat(determ,[1 size(Binv,2)]);

% Binv^T
BinvT = [Binv(:,1) Binv(:,3) Binv(:,2) Binv(:,4)];

determ = abs(determ);

% Compute quadrature nodes for all triangles on the physical domain
xx(1,:,:)=nodes_ref(:,1,:).*repmat(reshape(BE(:,1)',1,1,Ned),[Nqu,1,1])...
          +nodes_ref(:,2,:).*repmat(reshape(BE(:,3)',1,1,Ned),[Nqu,1,1])...
          +reshape(repmat(mesh.vertices1(E,1)',[Nqu 1]),Nqu,1,[]);
      
      
xx(2,:,:)=nodes_ref(:,1,:).*repmat(reshape(BE(:,2)',1,1,Ned),[Nqu,1,1])...
          +nodes_ref(:,2,:).*repmat(reshape(BE(:,4)',1,1,Ned),[Nqu,1,1])...
          +reshape(repmat(mesh.vertices1(E,2)',[Nqu 1]),Nqu,1,[]);

xx = permute(xx,[2 1 3]);


%Evaluate basis functions and their derivatives at quadrature points
vbasis=zeros(Nqu,Nloc,Ned);
dbasisx=zeros(Nqu,Nloc,Ned);
dbasisy=zeros(Nqu,Nloc,Ned);

switch Equation.base
    
    case 1
        
    % Monomial basis
    % ssn,ttn   : monomial values
    % dssn,dttn : monomial derivatives
    
    ssn=zeros(Nqu,degree+1,Ned);
    ttn=zeros(Nqu,degree+1,Ned);
    dssn=zeros(Nqu,degree+1,Ned);
    dttn=zeros(Nqu,degree+1,Ned);

    ssn(:,1,:)=1.0;
    ttn(:,1,:)=1.0;
    dssn(:,1,:)=0.0;
    dttn(:,1,:)=0.0;

    ssn(:,2,:)=nodes_ref(:,1,:);
    ttn(:,2,:)=nodes_ref(:,2,:);
    dssn(:,2,:)=1.0;
    dttn(:,2,:)=1.0;

    for i=3:degree+1
    
      ssn(:,i,:)=nodes_ref(:,1,:).*ssn(:,i-1,:);
      ttn(:,i,:)=nodes_ref(:,2,:).*ttn(:,i-1,:);
      dssn(:,i,:)=nodes_ref(:,1,:).*dssn(:,i-1,:);
      dttn(:,i,:)=nodes_ref(:,2,:).*dttn(:,i-1,:);
    
    end

    for i=3:degree+1

       dssn(:,i,:)= dssn(:,i,:)*(i-1);
       dttn(:,i,:)= dttn(:,i,:)*(i-1);
    
    end

    % values of basis functions 
    ii=1;
    for i=1:degree+1
       for j=1:i
        vbasis(:,ii,:)=ssn(:,i-j+1,:).*ttn(:,j,:);
        ii=ii+1;
       end
    end

   % derivatives of basis functions
   ii=1;
   for i=1:degree+1
       for j=1:i
        
        % local derivatives 
        dbasisx(:,ii,:)=dssn(:,i-j+1,:).*ttn(:,j,:);
        dbasisy(:,ii,:)=ssn(:,i-j+1,:).*dttn(:,j,:);

        ii=ii+1;
       end
   end

    case 2
        
    % Dubiner Basis
    
     ii=1;
        for i=1:degree+1
            for j=1:degree+1
                
                 if(i+j-2<=degree)
                 vbasis(:,ii,:)=Dubinerbasis(nodes_ref(:,1,:),nodes_ref(:,2,:),i-1,j-1);
                 [dbasisx(:,ii,:),dbasisy(:,ii,:)]=GradDubinerbasis(nodes_ref(:,1,:),nodes_ref(:,2,:),i-1,j-1);
                 ii=ii+1;
                 
                 end    
            end 
        end
        
        
      case  3
        % Only linear case 
        % 1-x-y, x, y are basis
        
       vbasis(:,1,:)=1-nodes_ref(:,1,:)-nodes_ref(:,2,:);
       vbasis(:,2,:)=nodes_ref(:,1,:);
       vbasis(:,3,:)=nodes_ref(:,2,:);
       
       dbasisx(:,1,:)=-ones(Nqu,1,Ned);
       dbasisx(:,2,:)=ones(Nqu,1,Ned);
       dbasisx(:,3,:)=zeros(Nqu,1,Ned);
       
       dbasisy(:,1,:)=-ones(Nqu,1,Ned);
       dbasisy(:,2,:)=zeros(Nqu,1,Ned);
       dbasisy(:,3,:)=ones(Nqu,1,Ned);
end

% Compute the derivatives dbasisx and dbasisy on physical element 
Pdbasisx=dbasisx.*repmat(reshape(BinvT(:,1)',1,1,Ned),[Nqu,Nloc,1])...
         +dbasisy.*repmat(reshape(BinvT(:,3)',1,1,Ned),[Nqu,Nloc,1]);

Pdbasisy=dbasisx.*repmat(reshape(BinvT(:,2)',1,1,Ned),[Nqu,Nloc,1])...
         +dbasisy.*repmat(reshape(BinvT(:,4)',1,1,Ned),[Nqu,Nloc,1]);

return;












  