%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%  By Murat UZUNCA and Hamdullah YÜCEL  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Main_Nonlinear()

% This routine solves the steady-state diffusion-convection-reaction equation
%
%   \alpha u - \epsilon*\Delta u + b\dot\nabla u + r(u) = f
%
% on the rectangular region [ax,bx]x[ay,by] using DG-FEM.

clear all
clc

%
%(ax,by)       (bx,by)
%   7-----8-----9 
%   |    /|    /|
%   |   / |   / |
%   |  /  |  /  |
%   | /   | /   |
%   4-----5-----6
%   |    /|    /|
%   |   / |   / |
%   |  /  |  /  |
%   | /   | /   |
%   1-----2-----3
%(ax,ay)       (bx,ay)

ax=0;    % x-coordinate of left boundary
bx=1;    % x-coordinate of right boundary
ay=0;    % y-coordinate of bottom boundary
by=1;    % y-coordinate of top boundary

% Generate the mesh
mx=0.5*(bx-ax);
my=0.5*(by-ay);    

% Nodes
Nodes = [ax,ay; mx,ay; bx,ay; ax,my; mx,my; bx,my; ax,by; mx,by;bx,by];                
% Elements (Counter-clock-wise direction)
Elements = [4,1,5; 1,2,5; 5,2,6; 2,3,6; 7,4,8; 4,5,8;8,5,9;5,6,9]; 
% Dirichlet bdry edges
Dirichlet = [1,2; 2,3; 1,4; 3,6; 4,7; 6,9; 7,8; 8,9];                
% Neumann bdry edges
Neumann   = []; 
%
% Initial mesh struct 
mesh = getmesh(Nodes,Elements,Dirichlet,Neumann);
dx=max(mx,my);
for jj=1:2
    mesh=uniformrefine(mesh);  %Refine mesh
    dx=0.5*dx;
end   

%% Setup Problem

% method : NIPG=1, SIPG=2, IIPG=3
method=2;
% Degree of polynomials
degree=1;  
    
% Set up the problem
[penalty,kappa]=set_parameter(method,degree);

% Existence of exact solution 
%0: No exact solution , 1: Exact solution exists
isexact=1;  

% Set up the number of levels and initialize the L2-error and dofs
% If exact solution exists
nlevel=4;
l2err=zeros(nlevel,1);
rate=zeros(nlevel,1);
dof=zeros(nlevel,1);
hmax=zeros(nlevel,1);

fprintf('    DoFs         Delta-x         L2-error         Rate         #it\n')

%% Solve Problem

% Loop over refinement levels
for level=1:nlevel

    %Compute global matrices and rhs global vector
    [D,C,R,F]=global_system(mesh,@fdiff,@fadv,@freact,...
        @fsource,@DBCexact,@NBCexact,penalty,kappa,degree);

    Stiff=D+C+R;   % Stiffness matrix

    % Initial guess for Newton iteration
    coef=zeros(size(Stiff,1),1);

    % Newton iteration
    noi=0;
    for ii=1:50
        noi=noi+1;

      % Compute the nonlinear vector and its Jacobian matrix at
      % the current iterate
       [H,HJ]=nonlinear_global(coef,mesh,@freact_nonlinear,degree);

      % Form the residual of the system
       Res = Stiff*coef + H - F;

      % Form the Jacobian matrix of the system 
      % (w.r.t. unknown coefficients coef)
       J = Stiff + HJ ; 

      % Solve the linear system for the correction "w"
       w = J \ (-Res);

      % Update the iterate
       coef = coef + w;

      % Check the accuracy
       if norm(J*w+Res) < 1e-20
           break;
       end

    end

    % Compute L2-error and plot the solution   
    if level==1
        [l2err(level),hmax(level)]=dg_error(coef,...
            mesh,@fexact,@fdiff,[],degree,isexact);
    else
        [l2err(level),hmax(level),rate(level)]=dg_error(coef,...
            mesh,@fexact,@fdiff,l2err(level-1),degree,isexact);
    end
    
    % Degree of freedom
    dof(level)=size(mesh.Elements,1)*(degree+1)*(degree+2)*0.5;    

    if level==1
        fprintf('%7d         %7s         %5.3e       %5s        %5d\n',...
            dof(level),strtrim(rats(dx)) ,l2err(level),'-',noi);
    else
        fprintf('%7d         %7s         %5.3e        %5.3f       %5d\n',...
            dof(level),strtrim(rats(dx)) ,l2err(level),rate(level),noi);
    end
    
    if level<nlevel
        mesh=uniformrefine(mesh);
        dx=0.5*dx;
    end
    
end

% Draw the errors
if ( nlevel>1 && isexact )
    draw_error(hmax,l2err,method,degree,nlevel,penalty)
end

end

%% Define diffusion, advection, and reaction as subfunctions

% Diffusion
function diff = fdiff(x,~)     
    diff = (10^(-3)).*ones(size(x));      
end

% Advection
function [adv1,adv2] = fadv(x,~) 
    adv1 =(1/sqrt(5))*ones(size(x)); 
    adv2 =(2/sqrt(5))*ones(size(x));
end

% Linear reaction 
function react = freact(x,~)
    react  = ones(size(x));
end

% Nonlinear reaction 
function [r,dr] = freact_nonlinear(u)
 % Value of the nonlinear reaction term at the current iterate
  r = u.^2;   

 % Value of the derivative of the nonlinear reaction 
 % term at the current iterate
  dr = 2*u;   
end

%% Define exact solution and force as subfunctions

% Exact solution
function [yex,yex_x,yex_y] = fexact(fdiff,x,y)    
    diff  = feval(fdiff,x,y);   % evaluate the diffusion function 
    yex=0.5*(1-tanh((2*x-y-0.25)./(sqrt(5*diff))));        % Exact value of adjoint variable
    yex_x=(-1./(sqrt(5*diff))).*(sech((2*x-y-0.25)./(sqrt(5*diff)))).^2;           % Derivative of adjoint variable wrt x.
    yex_y=((0.5)./(sqrt(5*diff))).*(sech((2*x-y-0.25)./(sqrt(5*diff)))).^2;    % Derivative of adjoint variable wrt y. 
end

% Force function
function source = fsource(fdiff,fadv,freact,x,y)    
    diff  = feval(fdiff,x,y );         % evaluate the diffusion function
    [adv1,adv2] = feval(fadv,x, y );   % evaluate the advection function 
    reac = feval(freact,x,y);        % evaluate the reaction function 

    yex=0.5*(1-tanh((2*x-y-0.25)./(sqrt(5*diff))));        % Exact value of adjoint variable
    yex_x=(-1./(sqrt(5*diff))).*(sech((2*x-y-0.25)./(sqrt(5*diff)))).^2;           % Derivative of adjoint variable wrt x.
    yex_y=((0.5)./(sqrt(5*diff))).*(sech((2*x-y-0.25)./(sqrt(5*diff)))).^2;    % Derivative of adjoint variable wrt y. 

    yex_xx=((0.8)./diff).*tanh((2*x-y-0.25)./(sqrt(5*diff))).*(sech((2*x-y-0.25)./(sqrt(5*diff)))).^2;
    yex_yy=((0.2)./diff).*tanh((2*x-y-0.25)./(sqrt(5*diff))).*(sech((2*x-y-0.25)./(sqrt(5*diff)))).^2;

    source=-diff.*(yex_xx+yex_yy)+(adv1.*yex_x+adv2.*yex_y)+reac.*yex+yex.^2;     
end


%% Boundary Conditions

% Drichlet Boundary Condition
function DBC=DBCexact(fdiff,x,y)
    % evaluate the diffusion function 
    diff  = feval(fdiff,x,y);   
    %Drichlet Boundary Condition
    DBC=0.5*(1-tanh((2*x-y-0.25)./(sqrt(5*diff))));                
end

% Neumann Boundary Condition
function NC = NBCexact(~,~,x,~)
    %Neumann Boundary Condition
    NC=zeros(size(x));      
end



%% Set-up  parameters function for DG-FEM

function [penalty,kappa]=set_parameter(method,degree)

global Equation;

Equation.b0=1;            % Superpenalization parameter (In standart b0=1, but to use superpenalization b0=3)
Equation.base=2;          % Choose the basis ( 1:monomials, 2:Dubiner Basis)

switch method
    case 1
        %NIPG
        Equation.method=1;   
        kappa=1;         % type of primal method
        penalty=1;     % penalty parameter
    case 2
        %SIPG
        Equation.method=2;
        kappa=-1;                        % type of primal method
        penalty=3*degree*(degree+1);   % penalty parameter
    case 3
        %IIPG
        Equation.method=3;
        kappa=0;                         % type of primal method
        penalty=3*degree*(degree+1);   % penalty parameter
   
end

end