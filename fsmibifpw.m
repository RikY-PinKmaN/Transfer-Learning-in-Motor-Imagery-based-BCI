function [fs,D,T]=fsmibifpw(d,varargin)
% FSMIBIFPW Feature Selection using MI based on Best Individual Features with Parzen Window
%   by Ang Kai Keng (kkang@pmail.ntu.edu.sg)
%
%   Syntax:
%     [fs,{D},{T}]=fsmibifpw(d,{nf})
%   where
%     fs : ranked indexes of top nf features.
%     D  : auxiliary results where
%          D.MC   = independent membership function classification results,
%          D.CO   = classification output,
%          D.SMC  = sorted membership classification results,
%          D.ISMC = sorted index of the membership classification results.
%     T  : time (s) taken to perform feature selection.
%     d  : data set matrix with size [sN,sD] where
%          sN          = number of data tuples,
%          sD-1        = number of input dimensions,
%          lastcol     = class id [0..nN-1].
%     nf : number of features to select where
%          (default) nf=sD-1.
%
%   See also MI, PARZENDE, FSCLASS.

% Determine the size of data matrix
[stN,stD]=size(d);
nf=stD-1;
pc.opt={{'nf', [], @isnumeric}};
pc.para={{'showwaitbar', false, @isnumeric},...
    {'nfgroup', 1, @isnumeric},...
    {'groups', 1, @isnumeric}};
parsearg(pc,varargin{:});

% Note start time
st=clock;

% Splits the data sets
otX=d(:,1:(stD-1));
otY=d(:,stD);

% Find the indexes to tX for each class
tXi0=find(otY==0);
tXi1=find(otY==1);

% Compute p(w_i)
pw0=length(tXi0)/stN;
pw1=length(tXi1)/stN;

% Initialize estimate p(x|w_i),p(w_i|x),H(w|x) and mi
pxw0=zeros(stN,stD-1);
pxw1=zeros(stN,stD-1);
pw0x=zeros(stN,stD-1);
pw1x=zeros(stN,stD-1);
hwx=zeros(1,stD-1);
D.mimf=zeros(1,stD-1);

% Compute Parzen Window width parameter
D.h=1/(log2(stN));

% Prepare the training data set
otX0=otX(tXi0,:);
otX1=otX(tXi1,:);

% Compute entropy of the class variable
hw=(-pw0*log2(pw0)-pw1*log2(pw1));

if showwaitbar
    h=waitbar(0,'FSMIBIFPW computing');
end

% Compute the parzen window classification results
for k=1:(stD-1)    
    if showwaitbar
        waitbar(k/(stD-1),h);
    end
    % Initialize conditional entropy H(w|x)
    hwx(k)=0;
    % Compute p(x_j|w_i) for each data point using parzen window
    pxw0(:,k)=parzende(otX0,k,otX(:,k));
    pxw1(:,k)=parzende(otX1,k,otX(:,k));
    t=pxw0(:,k)+pxw1(:,k);
    pw0x(:,k)=pxw0(:,k)./t;
    pw1x(:,k)=pxw1(:,k)./t;
    % Compute conditional entropy H(w|x) using all data tuples
    for j=1:stN
        % Compute estimate of conditional entropy H(w|x)
        if (pw0x(j,k)==0)
            hwx(k)=hwx(k)-(1/stN*(pw1x(j,k)*log2(pw1x(j,k))));
        else
            if (pw1x(j,k)==0)
                hwx(k)=hwx(k)-(1/stN*(pw0x(j,k)*log2(pw0x(j,k))));
            else
                hwx(k)=hwx(k)-(1/stN*(pw0x(j,k)*log2(pw0x(j,k))+pw1x(j,k)*log2(pw1x(j,k))));
            end
        end
    end;
    % Compute estimate of mi using feature k
    D.mimf(k)=hw-hwx(k);
end
if showwaitbar
    close(h);
end
% Group and sum the MI of each group
if groups==1
    mifg=vec2mat(1:stD-1,nfgroup)';
    mifsg=sum(vec2mat(D.mimf,nfgroup),2);
    % Sort the MI of each group
    [D.SMC D.ISMC]=sort(mifsg,'descend');
    % Return indexes of the top nf features
    fs=reshape(mifg(:,D.ISMC(1:nf)),1,[]);
else
    %mifg=vec2mat(1:stD-1,(stD-1)/groups);
    mifg=vec2mat(D.mimf,(stD-1)/groups);
    [D.SMC D.ISMC]=sort(mifg,2,'descend');
    [temp,D.segment]=max(sum(D.SMC(:,1:nf),2));
    fs={D.segment,D.ISMC(D.segment,1:nf)};
end    
% Compute time elapse
T=etime(clock,st);
