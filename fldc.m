function D=fldc(D,cd) 
% FLDC Fisher Linear Discriminant Compute function
%   by Ang Kai Keng (kkang@pmail.ntu.edu.sg)
%
%   Syntax:
%     D=fldc(D,cd) 
%   where
%     cd : compute data set matrix with size [sN,sD] where
%          sN          = number of data tuples,
%          sD-1        = number of input dimensions,
%          lastcol     = class id [0..nN-1].
%     D  : resultant data set where
%          D.pcO       = predicted computation output,
%          D.pcY       = predicted computation data class,
%          D.cconfuse  = computation confusion matrix,
%          D.caccuracy = computation classification accuracy.
%
%   See also FLDT.

% Local data
% ocX = original classification data,
% ocY = original classification data class,
% scX = standardized classification data,

% Start timer
st=clock;

% Determine the size of data matrix
[sN sD]=size(cd);

% Load X and Y matrix. D.oX=X', D.oY=Y'
ocX=cd(:,1:(sD-1));
ocY=cd(:,sD);

% Standardize data, D.csX=Standardized cX
scX = standardizew(ocX,D.net.sm,D.net.sd);

% compute mean of data and number of classes
m=mean(scX);
nN=1;

% Compute projection pX = z'*X.  
% Here D.pX=pX'=(z'*X)'=X'*z where scX=X', D.net.z=z
D.pcO=scX * D.net.z;

% Compute classification accuracy
D.pcY=zeros(sN,1);
if D.net.s==0
    for i=1:sN
        if (D.pcO(i,1)<D.net.b)
            D.pcY(i)=0;
        else
            D.pcY(i)=1;
        end
    end
else
    for i=1:sN
        if (D.pcO(i,1)>D.net.b)
            D.pcY(i)=0;
        else
            D.pcY(i)=1;
        end
    end   
end

% Compute confusion matrix
D.cconfuse=zeros(nN+1);
for i=1:nN+1;
    for j=1:nN+1;
        actual =(ocY==(i-1));
        predict=(D.pcY==(j-1));
        D.cconfuse(i,j)=sum(actual&predict);
    end;
end;

% Compute accuracy
D.caccuracy = sum(diag(D.cconfuse)/(sum(sum(D.cconfuse))))*100.0;

% Plot the projected data
% figure; dplot(D,1,2);
% 
% Plot discriminant line
% axis manual;
% plot([D.b D.b],[-9999 9999],'r');

% Stop timer
D.ct=etime(clock,st);
