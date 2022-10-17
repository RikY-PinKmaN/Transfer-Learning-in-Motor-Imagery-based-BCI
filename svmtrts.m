function [D,Dtest]=svmtrts(varargin)
% SVMT Support Vector Machine Train and Test function
%   by Ang Kai Keng (kkang@pmail.ntu.edu.sg)
%   Modified by Mahnaz Arvaneh (stuma@i2r.a-star.edu.sg)
%   Syntax:
%      [D,Dtest]=svmtrts(td,ts,'libsvm');
%   where
%     td : training data set [sN,sD] where
%          sN          = number of data tuples,
%          sD-1        = number of input dimensions,
%          lastcol     = class id [0..nN].
%     ts: Testing data set
%     D  : resultant data set where
%          D.net       = support vector machine
%          D.ptO       = predicted training output,
%          D.ptY       = predicted training data class,
%          D.tconfuse  = training confusion matrix,
%          D.taccuracy = training classification accuracy.
%
%   Parameter:
%     'standardize' : standardize training data before training
%                     (default) no standardize.
%     'showplot'    : show svmplot
%                     (default) svmplot is not shown.
%     'confi'       : compute distance D.pto from hyperplane
%                     (default) not computed.
%     'libsvm'      : use libsvm instead

% Local data
% otX = original training data,
% otY = original training data class.

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

%---------------------------------
%---------------------------------
if ~isempty(varargin)
    if isnumeric(varargin{1})
        test = varargin{1};
        varargin(1) = [];
    else
        error('Test data not specified!');
    end
else
    error('Test data not specified!');
end
%-------------------------------
%-------------------------------
showplot = false;
D.net.standardize = false;
confi = false;
libsvm = false;
% Check if parameters are specified
while ~isempty(varargin)
    switch varargin{1}
        case 'standardize'
            D.net.standardize = true;
            varargin(1) = [];
        case 'showplot'
            showplot = true;
            varargin(1) = [];
        case 'confi'
            confi = true;
            varargin(1) = [];
        case 'libsvm'
            libsvm = true;
            D.net.standardize = true;
            varargin(1) = [];
    end
end

% Determine the size of data matrix
[sN sD]=size(td);

% Splits the data set
otX=td(:,1:(sD-1));
otY=(td(:,sD));
%-------------------------------------
%-------------------------------------
[sNt sDt]=size(test);

% Splits the data set
otXt=test(:,1:(sDt-1));
otYt=(test(:,sDt));
%-------------------------------------
%-------------------------------------
% compute number of classes eg. 1=2 classes
nN=max(otY);

if D.net.standardize==true
    % Standardize data, D.osX=Standardized X
    [stX,D.net.sm,D.net.sd] = standardize(otX);
else
    stX=otX;
end
%-------------------------------------
%-------------------------------------
nNt=max(otYt);

if D.net.standardize==true
    % Standardize data, D.osX=Standardized X
    [stXt,D.net.sm,D.net.sd] = standardize(otXt);
else
    stXt=otXt;
end
%-------------------------------------
%-------------------------------------

% Train support vector machine
if libsvm==false
    options=optimset('MaxIter',1000);
    arg={'quadprog_opts',options};
    if showplot==true
        arg={arg{:},'showplot',true};
    end
    D.net.svm=svmtrain(stX,otY,arg{:});
else
    D.net.svm=libsvmtrain(otY,stX,'-t 0');
end
% Compute classification accuracy

% Evaluation training data ptY=class number starting from 1
if libsvm==false
    D.ptO=svmcompute(D.net.svm,stX);
    %D.ptO=svmclassify(D.net.svm,stX);
    D.ptY=(D.ptO<=0); % New version of svmclassify reversed the class labels
else
    D.ptO=libsvmpredict(otY,stX,D.net.svm);
    D.ptY=(D.ptO>=0);
end
%--------------------------
%--------------------------
%Evaluation testing data ptY=class number starting from 1
if libsvm==false
   Dtest.ptO=svmcompute(D.net.svm,stXt);
    %D.ptO=svmclassify(D.net.svm,stX);
    Dtest.ptY=(Dtest.ptO<=0); % New version of svmclassify reversed the class labels
else
    Dtest.ptO=libsvmpredict(otYt,stXt,D.net.svm);
    Dtest.ptY=(Dtest.ptO>=0);
end
%--------------------------
%--------------------------

% Compute confusion matrix
D.tconfuse=zeros(nN+1);
for i=1:nN+1;
    for j=1:nN+1;
        actual =(otY==(i-1));
        predict=(D.ptY==(j-1));
        D.tconfuse(i,j)=sum(actual&predict);
    end;
end;


%--------------------------
%--------------------------
% Compute confusion matrix
Dtest.tconfuse=zeros(nNt+1);
for i=1:nNt+1;
    for j=1:nNt+1;
        actualt =(otYt==(i-1));
        predictt=(Dtest.ptY==(j-1));
        Dtest.tconfuse(i,j)=sum(actualt&predictt);
    end;
end;

% Compute accuracy
D.taccuracy = sum(diag(D.tconfuse)/(sum(sum(D.tconfuse))))*100.0;

%----------------------------
%----------------------------
% Compute accuracy
Dtest.taccuracy = sum(diag(Dtest.tconfuse)/(sum(sum(Dtest.tconfuse))))*100.0;

% Stop timer
D.tt=etime(clock,st);
