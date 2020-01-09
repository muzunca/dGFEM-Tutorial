
% function [errL2,hmax,rate]=dg_error(coef,mesh,fexact,fdiff,err_ex,degree)
%
% Purpose : Compute the L2 error between the numerical solution 
%           and the exact function, and plot the solution 
% 
% Input   :
%          coef   : coefficients 
%          mesh   : mesh
%          fexact : exact function
%          fdiff  : diffusion function
%          err_ex : L2-error from previous level
%          degree : degree of polynomial functions
%
% Output  :
%          errL2  : L2 error between y and yex
%          hmax   : maximum h
%          rate   : convergence rate


function [errL2,hmax,rate]=dg_error(coef,mesh,fexact,fdiff,err_ex,degree,isexact)

set(0, 'defaultaxesfontsize',18,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',1.2,'defaultpatchlinewidth',1.2,...
       'defaulttextfontsize',18);
 
errL2=0;
rate=0;
   
Nloc=(degree+1)*(degree+2)/2;              % local dimension 
Nel=size(mesh.Elements,1);                 % number of elements

% Get quadrature points and weights on reference triangle
[nodes_ref,wgt]=quadrature(7);

% Compute the values of the basis functions,
% global derivative of the basis functions,
% the determinant of the transformation matrix and
% the values of quadrature points on all physical triangles
[val_basis,~,~,determ,xx]=elem_basis(mesh,degree,nodes_ref);

% quadrature weights
wgt=repmat(wgt,[1 1 Nel]);     % make Nel copies of the weights

ycoef =full(coef);                  % coefficients
yy    =reshape(ycoef,[Nloc,1,Nel]); % reshape coeff to Nloc x 1 x Nel
yval   = multiprod(val_basis,yy); 

if isexact
    % Evaluate yex and its gradient at the quadrature points
    [yex,~,~]=feval(fexact,fdiff,xx(:,1,:),xx(:,2,:));
    exvaly=squeeze(yex);     % Exact value

    % weights * determ and compute the transpose
    vol = permute(wgt.*determ,[2 1 3]);

    % L2 error
    errL2 = multiprod(vol, (yex - yval).^2);
    errL2 = sum(squeeze(errL2));
    errL2  = sqrt(errL2);

    if ~isempty(err_ex)
        rate=log2(err_ex/errL2)/log2(2);    
    end
end

% Compute hmax
adB    = determ(1,:,:);
adBmax = max(adB(:));
hmax   = (4*adBmax)^(1/2);

% Plot the exact and numerical solution
tri=delaunay(squeeze(xx(:,1,:)),squeeze(xx(:,2,:)));

numvaly=squeeze(yval);   % Numerical value

figure(700+degree)
if isexact
    subplot(1,2,1);
    trisurf(tri,squeeze(xx(:,1,:)),squeeze(xx(:,2,:)),exvaly(:))
    xlabel('x_{1}')
    ylabel('x_{2}')
    title('Exact Solution')
    shading 'interp';

    subplot(1,2,2);
end

trisurf(tri,squeeze(xx(:,1,:)),squeeze(xx(:,2,:)),numvaly(:))
shading 'interp';
xlabel('x_{1}')
ylabel('x_{2}')
title('DG Solution')

return;

end










