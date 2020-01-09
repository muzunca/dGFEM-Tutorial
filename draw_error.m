
% function draw_error(h,l2err,method,degree,nlevel,penalty)
%
% Purpose :  draw the graph of L2-error
%
% Input :   
%        h      : hmax (or DoFs)
%        l2err  : L2-error
%        method : DG method (NIPG=1,SIPG=2,IIPG=3)
%        degree : degree of polynomials
%        nlevel : number of refinement
%        penalty: penalty parameter


function draw_error(h,l2err,method,degree,nlevel,penalty)

global Equation

set(0, 'defaultaxesfontsize',14,'defaultaxeslinewidth',1.2,...
       'defaultlinelinewidth',1.2,'defaultpatchlinewidth',1.2,...
       'defaulttextfontsize',14);

% plot L2-error 
figure(method);
switch degree
    case 1
      loglog(h,l2err,'-rp', 'LineWidth',2); hold on 
      text(h(1)+0.1,l2err(1),'p=','fontsize',14), text(h(1)+0.21,l2err(1),num2str(degree),'fontsize',18)
      slope=log2(l2err(nlevel-1)/l2err(nlevel-0))/log2(2);
      plot([h(nlevel-0),h(nlevel-1),h(nlevel-1),h(nlevel-0)],...
           [l2err(nlevel-0),l2err(nlevel-1),l2err(nlevel-0),l2err(nlevel-0)],'r--')
      text(h(nlevel-1),(l2err(nlevel-0)+l2err(nlevel-1))*0.5,num2str(slope),'fontsize',18)

    case 2
      loglog(h,l2err,'--mo', 'LineWidth',2); hold on 
      text(h(1)+0.1,l2err(1),'p=','fontsize',18), text(h(1)+0.21,l2err(1),num2str(degree),'fontsize',18)
      slope=log2(l2err(nlevel-1)/l2err(nlevel-0))/log2(2);
      plot([h(nlevel-0),h(nlevel-1),h(nlevel-1),h(nlevel-0)],...
           [l2err(nlevel-0),l2err(nlevel-1),l2err(nlevel-0),l2err(nlevel-0)],'r--')
      text(h(nlevel-1),(l2err(nlevel-0)+l2err(nlevel-1))*0.5,num2str(slope),'fontsize',18)

    case 3
      loglog(h,l2err,':bs', 'LineWidth',2); hold on  
      text(h(1)+0.1,l2err(1),'p=','fontsize',18), text(h(1)+0.21,l2err(1),num2str(degree),'fontsize',18)
      slope=log2(l2err(nlevel-1)/l2err(nlevel-0))/log2(2);
      plot([h(nlevel-0),h(nlevel-1),h(nlevel-1),h(nlevel-0)],...
           [l2err(nlevel-0),l2err(nlevel-1),l2err(nlevel-0),l2err(nlevel-0)],'r--')
      text(h(nlevel-1),(l2err(nlevel-0)+l2err(nlevel-1))*0.5,num2str(slope),'fontsize',18)

    case 4
      loglog(h,l2err,'-.kd', 'LineWidth',2); hold on  
      text(h(1)+0.1,l2err(1),'p=','fontsize',18), text(h(1)+0.21,l2err(1),num2str(degree),'fontsize',18)
      slope=log2(l2err(nlevel-1)/l2err(nlevel-0))/log2(2);
      plot([h(nlevel-0),h(nlevel-1),h(nlevel-1),h(nlevel-0)],...
           [l2err(nlevel-0),l2err(nlevel-1),l2err(nlevel-0),l2err(nlevel-0)],'r--')
      text(h(nlevel-1),(l2err(nlevel-0)+l2err(nlevel-1))*0.5,num2str(slope),'fontsize',18)
      
    case 5
      loglog(h,l2err,'-g<', 'LineWidth',2); hold on  
      text(h(1)+0.1,l2err(1),'p=','fontsize',18), text(h(1)+0.21,l2err(1),num2str(degree),'fontsize',18)
      slope=log2(l2err(nlevel-1)/l2err(nlevel-0))/log2(2);
      plot([h(nlevel-0),h(nlevel-1),h(nlevel-1),h(nlevel-0)],...
           [l2err(nlevel-0),l2err(nlevel-1),l2err(nlevel-0),l2err(nlevel-0)],'r--')
      text(h(nlevel-1),(l2err(nlevel-0)+l2err(nlevel-1))*0.5,num2str(slope),'fontsize',18)
      
    case 6
      loglog(h,l2err,'--c>', 'LineWidth',2); hold on 
      text(h(1)+0.1,l2err(1),'p=','fontsize',18), text(h(1)+0.21,l2err(1),num2str(degree),'fontsize',18)
      slope=log2(l2err(nlevel-1)/l2err(nlevel-0))/log2(2);
      plot([h(nlevel-0),h(nlevel-1),h(nlevel-1),h(nlevel-0)],...
           [l2err(nlevel-0),l2err(nlevel-1),l2err(nlevel-0),l2err(nlevel-0)],'r--')
      text(h(nlevel-1),(l2err(nlevel-0)+l2err(nlevel-1))*0.5,num2str(slope),'fontsize',18)
      
    case 7
      loglog(h,l2err,':yh', 'LineWidth',2); hold on  
      text(h(1)+0.1,l2err(1),'p=','fontsize',18), text(h(1)+0.21,l2err(1),num2str(degree),'fontsize',18)
      slope=log2(l2err(nlevel-1)/l2err(nlevel-0))/log2(2);
      plot([h(nlevel-0),h(nlevel-1),h(nlevel-1),h(nlevel-0)],...
           [l2err(nlevel-0),l2err(nlevel-1),l2err(nlevel-0),l2err(nlevel-0)],'r--')
      text(h(nlevel-1),(l2err(nlevel-0)+l2err(nlevel-1))*0.5,num2str(slope),'fontsize',18)
        
end

xlabel('DoFs');
ylabel('L_{2} Error');  
switch method
     case 1
         if(Equation.b0~=1)
         title('NIPG with superpenalization ','fontsize',20);
         elseif (penalty==0)
         title('NIPG0 Method ','fontsize',20);   
         else
         title('NIPG Method ','fontsize',20); 
         end
     case 2
         title('SIPG Method','fontsize',20);
     case 3
         if(Equation.b0~=1)
         title('IIPG with superpenalization','fontsize',20);
         else
         title('IIPG Method','fontsize',20);   
         end
end

return;















