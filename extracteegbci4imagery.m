function [data,temp]=extracteegbci4imagery(raweeg,xmethod,varargin)
% EXTRACTEEGDATA Extract data from Raw EEG data
%   by Ang Kai Keng (kkang@i2r.a-star.edu.sg)
%
%   Syntax:
%     data=extracteegdata(raweeg)
%     data=extracteegdata(raweeg,xmethod)
%   where
%     data:    EEG Data structure for classification
%     raweeg:  Raw EEG data
%     xmethod: 'depressed','pressed','indicate','all'
%     'depressed' extracts EEG data at key-depressed event
%     'pressed'   extracts EEG data at key-pressed event
%     'indicate'  extracts EEG data at indication event
%     'all'       extracts No event and key-pressed event
%     'specific'  extracts specified stimcode
%
%   Parameter:
%     seconds  : specify the number of seconds of data before and after
%                event to extract (default) [-2 2]
%     npress   : specify which key-pressed event to extract (default) 1-3
%     selchs   : specify which channels to extract (default 15 channels)
%     relax    : extract relax between action events
%     stim     : specify stimcode to extract
%
%   See also PREPEEGDATA, LOADEEGDATA.

% Default 'pressed' extraction method
% if nargin<2
%     xmethod='pressed';
% end
% Default nsec to extract is 2
nsec=[-1 1];
if strcmpi(xmethod,'all')
    nsec=[0 4];
end
% Default extracts all 3 keypress
npress=0; 
extractrelax=0;

% Default extracts only 15 channels from all 40
%strChList={'F3','Fz','F4','FC3','FCz','FC4','C3','Cz','C4','CP3','CPz','CP4','P3','Pz','P4'};
strChList=raweeg.sel_chan_list;
eeglen=size(raweeg.EEG,2);
detectoutlier=false;
% Process parameters
while ~isempty(varargin)
    if ischar(varargin{1})
        switch varargin{1}
            case 'npress'
                npress=varargin{2};
                varargin(1)=[];
                varargin(1)=[];            
            case 'seconds'
                nsec=varargin{2};
                varargin(1)=[];
                varargin(1)=[];
            case 'selchs'
                strChList=varargin{2};
                varargin(1)=[];
                varargin(1)=[];
            case 'relax'
                extractrelax=1;
                varargin(1)=[];
            case 'relax2'
                extractrelax=2;
                varargin(1)=[];
            case 'relax3'
                extractrelax=3;
                varargin(1)=[];
            case 'stim'
                specstim=varargin{2};
                varargin(1)=[];
                varargin(1)=[];                
            otherwise
                break;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCh=length(strChList);
IndChs=zeros(1,nCh);
for iCh=1:nCh
    for id=1:size(raweeg.EEG,1)
        if(strcmpi(strChList{iCh},raweeg.chan_list{id})) 
            IndChs(iCh)=id;
            break
        end
    end
end

% Identify outliers, which is defined as deviations from the range [-5000,5000]
if detectoutlier
    outlier_raweeg=abs(raweeg.EEG)>5000;
else
    outlier_raweeg=isnan(raweeg.EEG);
end

% Interpretation for the stimcode from raweeg data without deleting bad trials
start_trial_pos1    = raweeg.trial_start_pos;
end_trial_pos1       = raweeg.trial_end_pos;
trial_cue_pos1       =raweeg.trial_cue_pos; 
% cue
indicate_left_pos1   = raweeg.stimpos((raweeg.stimcode==769));
indicate_right_pos1  = raweeg.stimpos((raweeg.stimcode==770));
trial_class1= raweeg.trial_class;
trial_bad_mark=raweeg.trial_bad_mark;

% Interpretation for the stimcode from raweeg data after deleting bad trials

start_trial_pos    = raweeg.trial_start_pos(trial_bad_mark==0);
end_trial_pos       = raweeg.trial_end_pos(trial_bad_mark==0);
trial_cue_pos       =raweeg.trial_cue_pos(trial_bad_mark==0);  
%cue
indicate_left_pos   = raweeg.trial_cue_pos((raweeg.trial_class==1& trial_bad_mark==0));
indicate_right_pos  = raweeg.trial_cue_pos((raweeg.trial_class==2& trial_bad_mark==0));
trial_class= raweeg.trial_class(trial_bad_mark==0);
trial_bad_mark=raweeg.trial_bad_mark;

% Compute number of data points to extract
ptsrange=nsec*raweeg.sampling_rate;
ptsperCh=(ptsrange(2)-ptsrange(1))+1;
        

%===========================================================
 if ~strcmp(xmethod,'specific')
    % Initialize local variables
    nleft=0;
    nright=0;
    nleft=0;
    nright=0;
    left_start_pos     = zeros(size(indicate_left_pos));
    left_indicate_pos  = zeros(size(indicate_left_pos));
    %left_cue_pos   = zeros(3,length(trial_cue_pos(raweeg.trial_class==769)));
    left_cue_pos   = zeros(3,length(indicate_left_pos));
    left_end_pos       = zeros(size(indicate_left_pos));
    right_start_pos    = zeros(size(indicate_right_pos));
    right_indicate_pos = zeros(size(indicate_right_pos));
    right_pressed_pos  = zeros(3,length(right_start_pos));
    right_end_pos      = zeros(size(indicate_right_pos));
    left_outlier       = false(size(indicate_left_pos));
    right_outlier      = false(size(indicate_right_pos));

    % Distinguish between a trial for the left or right
%     for trial=1:length(start_trial_pos)
%        indicate_left=find(indicate_left_pos>start_trial_pos(trial) & indicate_left_pos<end_trial_pos(trial));
%         indicate_right=find(indicate_right_pos>start_trial_pos(trial) & indicate_right_pos<end_trial_pos(trial));
%      
%        
%         if ~isempty(indicate_left)
%             nleft=nleft+1;
%             left_start_pos(nleft)=indicate_left_pos(trial);
%            %left_end_pos(nleft)=end_trial_pos(trial);
%             left_indicate_pos(nleft)=indicate_left_pos(indicate_left);
% %             if ~strcmpi(xmethod,'depressed') && detectoutlier
% %                 if length(pressed_left)<3
% %                     disp(['Outlier: Less than 3 key pressed in ' num2str(nleft) ' left trial']);
% %                     left_outlier(nleft)=true;
% %                     if length(pressed_left)<npress
% %                         error(['There is no ' num2str(npress) ' in data']);
% %                     else
% %                         left_pressed_pos(:,nleft)=[pressed_left_pos(pressed_left) zeros(1,3-length(pressed_left))];
% %                     end
% %                 else
% %                     left_pressed_pos(:,nleft)=pressed_left_pos(pressed_left);                
% %                 end
% %             end
% %             left_end_pos(nleft)=end_trial_pos(trial);
%         else
%             if ~isempty(indicate_right)
%                 nright=nright+1;
%                 right_start_pos(nright)=start_trial_pos(trial);
%                 right_end_pos(nright)=end_trial_pos(trial);
%                 right_indicate_pos(nright)=indicate_right_pos(indicate_right);
%                 if ~strcmpi(xmethod,'depressed') && detectoutlier
%                     if length(pressed_right)<3
%                         disp(['Outlier: Less than 3 key pressed in ' num2str(nright) ' right trial']);
%                         right_outlier(nright)=true;
%                         if length(pressed_right)<npress
%                             error(['There is no ' num2str(npress) ' in data']);
%                         else
%                             right_pressed_pos(:,nright)=[pressed_right_pos(pressed_right) zeros(1,3-length(pressed_right))];
%                         end
%                     else
%                         right_pressed_pos(:,nright)=pressed_right_pos(pressed_right);                    
%                     end 
%                 end
%                 right_end_pos(nright)=end_trial_pos(trial);
%             else
%                 error('No indication found between start and end of trial');
%             end
%         end    
%     end
%==========================================================================
 end
% case 'indicate'
        %data.EEG_Left{1}=zeros(nCh,ptsperCh,nleft);
        %data.EEG_Right{1}=zeros(nCh,ptsperCh,nright);
        nleft=length(indicate_left_pos);
        nright=length(indicate_right_pos);
        EData.Left=zeros(ptsperCh,nCh,nleft);
        EData.Right=zeros(ptsperCh,nCh,nright);
        EData.x=zeros(ptsperCh,nCh,nleft+nright);
        EData.y=zeros(nleft+nright,1);
        EData.Relax=[];
        for i=1:nleft
            data_idx=indicate_left_pos(i)+ptsrange(1):indicate_left_pos(i)+ptsrange(2);
%             if left_outlier(i)
%                 data_left_outlier(i)=true;
%                 disp(['Skip: Outlier in ' num2str(i) ' of left indicate']);                        
%             end                    
            for j=1:nCh
                EData.Left(:,j,i)=[zeros(1,length(data_idx(data_idx<=0)))...
                                   double(raweeg.EEG(IndChs(j),data_idx(data_idx>0 & data_idx<=eeglen)))*raweeg.resolution...
                                   zeros(1,length(data_idx(data_idx>eeglen)))];
            end
        end
        for i=1:nright
            data_idx=indicate_right_pos(i)+ptsrange(1):indicate_right_pos(i)+ptsrange(2);
%             if right_outlier(i)
%                 data_right_outlier(i)=true;
%                 disp(['Skip: Outlier in ' num2str(i) ' of right indicate']);
%             end
            for j=1:nCh
                EData.Right(:,j,i)=[zeros(1,length(data_idx(data_idx<=0)))...
                                    double(raweeg.EEG(IndChs(j),data_idx(data_idx>0 & data_idx<=eeglen)))*raweeg.resolution...
                                    zeros(1,length(data_idx(data_idx>eeglen)))];
            end
        end
        EData.x(:,:,1:nleft)= EData.Left;
        EData.x(:,:,1+nleft:nright+nleft)=EData.Right;
        EData.y(1:size(EData.Left,3),1)=1;
        EData.y(size(EData.Left,3)+1:size(EData.Right,3)+size(EData.Left,3),1)=2;
        %sampling rate
        EData.s=250;
        data=EData;
%         if nrelax>0
%             EData.Relax=zeros(ptsperCh,nCh,nrelax);
%             for i=1:nrelax
%                 data_idx=indicate_relax_pos(i)+ptsrange(1):indicate_relax_pos(i)+ptsrange(2);
%                 % No outlier check for relaxed
%                 for j=1:nCh
%                     %data.EEG_Right{1}(j,:,i)=raweeg.EEG(IndChs(j),data_idx);
%                     EData.Relax(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
%                 end                
%             end
%         end
