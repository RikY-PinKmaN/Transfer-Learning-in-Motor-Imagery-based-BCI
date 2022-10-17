function [W Cvar]=cspbci4(xdata,varargin)
% CSP Computes the CSP Projection Matrix W from EEG data
%   by Ang Kai Keng (kkang@i2r.a-star.edu.sg)
%
%   Syntax:
%     W=csp(xdata)
%   where
%     xdata: extracted eegdata
%
%   See also extracteegdata.

nch=size(xdata.x,2);
ntrials=size(xdata.x,3);
nclass=2;
leavesingle=false;
eigmethod = 'eig';
covmethod = 'Ramoser';
compute_var=false;

if nargout==2
    compute_var=true;
    Cvar=[];
end
while ~isempty(varargin)
    if ischar(varargin{1})
        switch varargin{1}
            case 'leavesingle'
                leavesingle=true;
                varargin(1)=[];
            case 'svd'
                eigmethod='svd';
                varargin(1)=[];
            case 'eig'
                eigmethod='eig';
                varargin(1)=[];
            case 'blankertz';
                covmethod='Blankertz';
                varargin(1)=[];
            otherwise
                varargin(1)=[];
        end
    end
end

% Prealloc variables
C=zeros(nch,nch,ntrials);
Cmean=cell(nclass);
for i=1:nclass
    Cmean{i}=zeros(nch,nch);
end

switch covmethod
    case 'Ramoser'
        % Ramoser equation (1)
        for trial=1:ntrials
            E=xdata.x(:,:,trial)'; % Transpose samples x channels to channels x samples
            %E=E-repmat(mean(E,2),[1 nsamples]); % Zero mean
            tmpC = (E*E');
            C(:,:,trial) = tmpC./trace(tmpC);
            %C(:,:,trial) = tmpC;
        end
        for i=1:nclass
            Cmean{i}=mean(C(:,:,(xdata.y==i)),3);
            if compute_var
                Cvar{i}=var(C(:,:,(xdata.y==i)),0,3);
            end
        end
    case 'Blankertz'
        for i=1:nclass
            E=num2cell(xdata.x(:,:,xdata.y==(i)),[1 2]);
            X=cat(1,E{:});
            tmpC=X'*X;
            Cmean{i}=tmpC./trace(tmpC);
        end
end

% Ramoser equation (2)
Ccompo=Cmean{1}+Cmean{2};

% Find rank of Ccompo (new step to eliminate matrix singularity at line 53)
if leavesingle
    Crank=max(size(Ccompo));
else
    Crank=rank(Ccompo);
end


switch eigmethod
    case 'eig'
        [Ualter,Lalter]=eig(Ccompo,Cmean{1});
        % Sort eigenvalues and eigenvectors
        [Lalter,ind]=sort(diag(Lalter),'descend');
        Ualter=Ualter(:,ind);
        % Retain up to rank eigenvalues and eigenvectors (new step to eliminate matrix singularity at line 53)
        Ualter=Ualter(:,1:Crank);
        W=Ualter';
    case 'svd'
        [U,S,V] = svd(Ccompo);    %#ok<NASGU> %the eigenvalues are in decreasing order
        eigenvalues = diag(S);
        eigenvalues = eigenvalues(1:Crank);
        U = U(:,1:Crank);
        eigenvalues = sqrt(eigenvalues);
        eigenvalues = 1./eigenvalues;
        P=diag(eigenvalues)*U';
        S1 = P*Cmean{1}*P';
        %S2 = P*Cmean{2}*P';
        [u1,eigen1,v1]=svd(S1); %#ok<NASGU>
        %[ub,eigenB,vb]=svd(S2);
        W = u1'*P; % get the projection matrix
end

%Normalize projection matrix W
for i=1:length(W) %W is always a square matrix
    W(i,:)=W(i,:)./norm(W(i,:));
end
% Common spatial patterns in inverse of W
%A=pinv(W); % Common spatial patterns
