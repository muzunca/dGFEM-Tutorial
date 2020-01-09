

% function Elements = label(Nodes,Elements)
%
% Purpose : label the longest edge of each triangle as the base
% 
% Input:
%       Nodes    :  coordinates of Nodess
%       Elements :  Elementsent array
%
% Output:
%       Elements :  Elements(t,1) is opposite to the longest edge
%
%  This routine is a modified form of AFEM by L. Chen & C. Zhang 


function Elements = label(Nodes,Elements)


%--------------------------------------------------------------------------
% Compute length of each edge
%--------------------------------------------------------------------------
edgelength = zeros(size(Elements)); % initialize to accelerate access

edgelength(:,1) = (Nodes(Elements(:,3),1)-Nodes(Elements(:,2),1)).^2 ...
                + (Nodes(Elements(:,3),2)-Nodes(Elements(:,2),2)).^2;
edgelength(:,2) = (Nodes(Elements(:,1),1)-Nodes(Elements(:,3),1)).^2 ...
                + (Nodes(Elements(:,1),2)-Nodes(Elements(:,3),2)).^2;
edgelength(:,3) = (Nodes(Elements(:,1),1)-Nodes(Elements(:,2),1)).^2 ...
                + (Nodes(Elements(:,1),2)-Nodes(Elements(:,2),2)).^2;

%--------------------------------------------------------------------------
% Re-labelling according the longest edge
%--------------------------------------------------------------------------
[~,I] = max(edgelength,[],2);
Elements((I==2),[1 2 3]) = Elements((I==2), [2 3 1]);
Elements((I==3),[1 2 3]) = Elements((I==3), [3 1 2]);

return;



