% RCSP optimaization function
%   by Mahnaz Arvaneh(stuma@i2r.a-star.edu.sg)
function f = myfun(x)
global sigmac
global sigmat
global C
d=22;
% sumrayl=0.0001;
% for i=1:22
% rayl= (x(:,i)'*(sigmac)*x(:,i))/(x(:,i)'*sigmat*x(:,i));
% modifterm=C*norm(x(:,i),1)/((d^0.5)*norm(x(:,i)));
% %modifterm=C*norm(x(:,i),2)/((d^0.5)*norm(x(:,i),1));
% if (rayl-modifterm)>0
% sumrayl=sumrayl+1/(rayl-modifterm);
% else
%     sumrayl=sumrayl;
% end
%end
f=x*sigmac*x'+C*norm(x,1)/norm(x,2);
