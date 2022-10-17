function D=fldt(varargin) 
% FLDT Fisher Linear Discriminant Train function
%   by Ang Kai Keng (kkang@pmail.ntu.edu.sg)
%
%   Syntax:
%     D=fldt(td) 
%   where
%     td : training data set matrix with size [sN,sD] where
%          sN          = number of data tuples,
%          sD-1        = number of input dimensions,
%          lastcol     = class id [0..nN-1].
%     D  : resultant data set where
%          D.net.sm    = mean of training data,
%          D.net.sd    = standard deviation of training data,
%          D.net.z     = transformation matrix,
%          D.net.b     = discriminant bias,
%          D.ptO       = predicted training output,
%          D.ptY       = predicted training data class,
%          D.tconfuse  = training confusion matrix,
%          D.taccuracy = training classification accuracy.
%
%   See also FLDC.

% Local data
% otX = original training data,
% otY = original training data class,
% osX = standardized training data,

% Start timer
st=clock;

% Check if D is specified
if ~isempty(varargin)
    if isstruct(varargin{1}) || isempty(varargin{1})
        D = varargin{1};
        varargin(1) = [];
    end
end

% Check if td is specified
if ~isempty(varargin)
    if isnumeric(varargin{1})
        td = varargin{1};
        varargin(1) = [];
    else
        error('Training data not specified!');
    end
else
    error('Training data not specified!');
end

% Determine the size of data matrix
[sN sD]=size(td);

% Load X and Y matrix. D.oX=X', D.oY=Y'
otX=td(:,1:(sD-1));
otY=td(:,sD);

% Standardize data, D.osX=Standardized X
[osX,D.net.sm,D.net.sd] = standardize(otX);

% compute mean of data and number of classes
m=mean(osX);
nN=max(otY);
if nN>1
    error('Supports only 2 class computation!');
end

% initialize Sw (within class scatter) Sb (between class scatter)
Sw=zeros(sD-1,sD-1);
Sb=zeros(sD-1,sD-1);

% Compute Sb and Sw
for i=0:nN
	% split data into each class
	xi=osX(find(otY==i),:);
	% compute each mean (row vector)
	mi=mean(xi);
	% compute each cov
    Syi=cov(xi)*(length(xi)-1);
	% compute Sb and Sw. Here (mi-m)' is (u_yi-u)
    Sw=Sw+Syi;
	Sb=Sb+length(xi)*(mi-m)'*(mi-m);
end;

% Compute mean for bias computation
x0=osX(find(otY==0),:);
x1=osX(find(otY==1),:);
m0=mean(x0);
m1=mean(x1);

% Compute D.C=inv(Sw)*Sb
C=pinv(Sw)*Sb;

% Visualize first 2 components
[D.net.z,l,e] = pcacov(C);

% Formula for projection pX = z'*X.  
% Here D.ptO=pX'=(z'*X)'=X'*z where osX=X', D.net.z=z
D.ptO  = osX * D.net.z;

% Compute bias
m0=(m0*D.net.z);
m1=(m1*D.net.z);
D.net.b = (m1(1)+m0(1))/2;
D.net.s = 0;

% Compute classification accuracy
D.ptY=zeros(sN,1);
for i=1:sN
    if (D.ptO(i,1)<D.net.b)
        D.ptY(i)=0;
    else
        D.ptY(i)=1;
    end;
end;
if (length(find(D.ptY==otY))<(0.5*sN))
    D.net.s = 1;
    pY0=find(D.ptY==0);
    pY1=find(D.ptY==1);
    D.ptY(pY0)=1;
    D.ptY(pY1)=0;
end

% Compute confusion matrix
D.cconfuse=zeros(nN+1);
for i=1:nN+1;
    for j=1:nN+1;
        actual =(otY==(i-1));
        predict=(D.ptY==(j-1));
        D.tconfuse(i,j)=sum(actual&predict);
    end;
end;

% Compute accuracy
D.taccuracy = sum(diag(D.tconfuse)/(sum(sum(D.tconfuse))))*100.0;

% Stop timer
D.tt=etime(clock,st);

% Plot the projected data
% figure; dplot(D,1,2);
%
% Plot discriminant line
% axis manual;
% plot([D.b D.b],[-9999 9999],'r');