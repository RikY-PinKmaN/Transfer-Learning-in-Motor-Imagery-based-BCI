function [data,temp]=extracteegdata(raweeg,xmethod,varargin)
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
if nargin<2
    xmethod='pressed';
end
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
strChList = [];
if isfield(raweeg, 'sel_chan_list')
    strChList = raweeg.sel_chan_list;
end

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
            case 'stim'
                specstim=varargin{2};
                varargin(1)=[];
                varargin(1)=[];                
            otherwise
                break;
        end
    end
end

if isempty(strChList)
    if isfield(raweeg, 'ch_idx')
        IndChs = raweeg.ch_idx + 1;
        nCh = length(IndChs);
    end
else
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
end

if isempty(IndChs)
    error('Used channels not defined!');
end

% Identify outliers, which is defined as deviations from the range [-5000,5000]
if detectoutlier
    outlier_raweeg=abs(raweeg.EEG)>5000;
else
    outlier_raweeg=isnan(raweeg.EEG);
end

% Interpretation for the stimcode from raweeg data
start_trial_pos     = raweeg.stimpos((raweeg.stimcode==100));
end_trial_pos       = raweeg.stimpos((raweeg.stimcode==199));
indicate_relax_pos  = raweeg.stimpos((raweeg.stimcode==120));    
indicate_left_pos   = raweeg.stimpos((raweeg.stimcode==121));
indicate_right_pos  = raweeg.stimpos((raweeg.stimcode==122));
depressed_left_pos  = raweeg.stimpos((raweeg.stimcode==131));
depressed_right_pos = raweeg.stimpos((raweeg.stimcode==132));
pressed_left_pos    = raweeg.stimpos((raweeg.stimcode==151));
pressed_right_pos   = raweeg.stimpos((raweeg.stimcode==152));

% Compute number of data points to extract
ptsrange=floor(nsec*raweeg.sampling_rate);
ptsperCh=(ptsrange(2)-ptsrange(1))+1;

if ~strcmp(xmethod,'specific')
    % Initialize local variables
    nleft=0;
    nright=0;
    left_start_pos     = zeros(size(indicate_left_pos));
    left_indicate_pos  = zeros(size(indicate_left_pos));
    left_pressed_pos   = zeros(3,length(left_start_pos));
    left_end_pos       = zeros(size(indicate_left_pos));
    right_start_pos    = zeros(size(indicate_right_pos));
    right_indicate_pos = zeros(size(indicate_right_pos));
    right_pressed_pos  = zeros(3,length(right_start_pos));
    right_end_pos      = zeros(size(indicate_right_pos));
    left_outlier       = false(size(indicate_left_pos));
    right_outlier      = false(size(indicate_right_pos));

    % Distinguish between a trial for the left or right
    for trial=1:length(start_trial_pos)
        indicate_left=find(indicate_left_pos>start_trial_pos(trial) & indicate_left_pos<end_trial_pos(trial));
        indicate_right=find(indicate_right_pos>start_trial_pos(trial) & indicate_right_pos<end_trial_pos(trial));
        pressed_left=find(pressed_left_pos>start_trial_pos(trial) & pressed_left_pos<end_trial_pos(trial));
        pressed_right=find(pressed_right_pos>start_trial_pos(trial) & pressed_right_pos<end_trial_pos(trial));   
        if ~isempty(indicate_left)
            nleft=nleft+1;
            left_start_pos(nleft)=start_trial_pos(trial);
            left_end_pos(nleft)=end_trial_pos(trial);
            left_indicate_pos(nleft)=indicate_left_pos(indicate_left);
            if ~strcmpi(xmethod,'depressed') && detectoutlier
                if length(pressed_left)<3
                    disp(['Outlier: Less than 3 key pressed in ' num2str(nleft) ' left trial']);
                    left_outlier(nleft)=true;
                    if length(pressed_left)<npress
                        error(['There is no ' num2str(npress) ' in data']);
                    else
                        left_pressed_pos(:,nleft)=[pressed_left_pos(pressed_left) zeros(1,3-length(pressed_left))];
                    end
                else
                    left_pressed_pos(:,nleft)=pressed_left_pos(pressed_left);                
                end
            end
            left_end_pos(nleft)=end_trial_pos(trial);
        else
            if ~isempty(indicate_right)
                nright=nright+1;
                right_start_pos(nright)=start_trial_pos(trial);
                right_end_pos(nright)=end_trial_pos(trial);
                right_indicate_pos(nright)=indicate_right_pos(indicate_right);
                if ~strcmpi(xmethod,'depressed') && detectoutlier
                    if length(pressed_right)<3
                        disp(['Outlier: Less than 3 key pressed in ' num2str(nright) ' right trial']);
                        right_outlier(nright)=true;
                        if length(pressed_right)<npress
                            error(['There is no ' num2str(npress) ' in data']);
                        else
                            right_pressed_pos(:,nright)=[pressed_right_pos(pressed_right) zeros(1,3-length(pressed_right))];
                        end
                    else
                        right_pressed_pos(:,nright)=pressed_right_pos(pressed_right);                    
                    end 
                end
                right_end_pos(nright)=end_trial_pos(trial);
            else
                error('No indication found between start and end of trial');
            end
        end    
    end

    % Compute relax position
    if isempty(indicate_relax_pos) && extractrelax>0
        switch xmethod
            case 'pressed'
                error('Feature not implemented yet');
            case 'indicate'
                indicate_pos  = raweeg.stimpos((raweeg.stimcode==121) | (raweeg.stimcode==122));
                indicate_relax_pos=zeros(1,length(indicate_pos));
                indicate_pos = [indicate_pos length(raweeg.EEG)]; % Add last pos
                len=length(indicate_relax_pos);
                switch extractrelax
                    case 1
                        for i=1:len
                            bound=max(abs(nsec))*raweeg.sampling_rate;
                            relax_len=length(ptsrange(1):ptsrange(2));
                            start_relax_pos=indicate_pos(i)+bound;
                            end_relax_pos=indicate_pos(i+1)-bound;
                            if start_relax_pos<(end_relax_pos-relax_len)
                                indicate_relax_pos(i)=start_relax_pos+floor((end_relax_pos-start_relax_pos-relax_len)/2);
                            end
                        end
                    case 2
                        % Specific for BCI Competition III dataset IVc

                        event_interval=5.5*raweeg.sampling_rate;

                        % Estimate number of relax (tongue) events between indicate intervals
                        nrelax_interval=round(((indicate_pos(2:len+1)-indicate_pos(1:len)))/event_interval)-1;

                        % Delete number of relax (tongue) envents if it is more than 4
                        nrelax_interval(nrelax_interval>4)=0;

                        relax_len=sum(nrelax_interval);
                        indicate_relax_pos=zeros(1,relax_len);
                        relax_pos_idx=0;
                        for i=1:len
                            for j=1:nrelax_interval(i)
                                relax_pos_idx=relax_pos_idx+1;
                                indicate_relax_pos(relax_pos_idx)=j*event_interval+indicate_pos(i);
                            end
                        end
                end
                indicate_relax_pos(indicate_relax_pos==0)=[];
        end
        %temp=nrelax_interval;
        temp=indicate_relax_pos;
    end

    % Detect outlier trial
    for i=1:nleft
        data_idx=left_start_pos(i):left_end_pos(i);
        if sum(sum(outlier_raweeg(IndChs(:),data_idx)))>1
            left_outlier(i)=true;
        end
    end
    for i=1:nright
        data_idx=right_start_pos(i):right_end_pos(i);
        if sum(sum(outlier_raweeg(IndChs(:),data_idx)))>1
            right_outlier(i)=true;
        end
    end   
    clear nleft;
    clear nright;

    switch xmethod
        case 'depressed'
            % Error checking        
            nleft=length(indicate_left_pos);
            nright=length(indicate_right_pos);

            if length(depressed_left_pos)>nleft
                disp('Warning: Extra left pressed detected!');
                i=1;
                while i<nleft
                    if depressed_left_pos(i)<left_start_pos(i) || depressed_left_pos(i)>left_end_pos(i)                  
                        depressed_left_pos(i)=[];
                        disp('Warning: Extra left pressed removed.');
                    else
                        i=i+1;
                    end
                end            
            end
            if length(depressed_right_pos)>nright
                disp('Warning: Extra right pressed detected!');
                i=1;
                while i<nright
                    if depressed_right_pos(i)<right_start_pos(i) || depressed_right_pos(i)>right_end_pos(i)                  
                        depressed_right_pos(i)=[];
                        disp('Warning: Extra right pressed removed.');
                    else
                        i=i+1;
                    end
                end            
            end

        case 'pressed'       
            if npress==0
                nleft=length(pressed_left_pos);
                nright=length(pressed_right_pos);
            else
                nleft=length(indicate_left_pos);
                nright=length(indicate_right_pos);        
            end
        case {'indicate','all'}
            nleft=length(indicate_left_pos);
            nright=length(indicate_right_pos);
            nrelax=length(indicate_relax_pos);
        otherwise
            nleft=0;
            nright=0;
    end        

    data_left_outlier=false(1,nleft);
    data_right_outlier=false(1,nright);
end

switch xmethod
    case 'depressed'
        EData.Left =zeros(ptsperCh,nCh,nleft);
        EData.Right=zeros(ptsperCh,nCh,nright);
        EData.Relax=[];
        for i=1:nleft
            data_idx=depressed_left_pos(i)+ptsrange(1):depressed_left_pos(i)+ptsrange(2);
            for j=1:nCh                    
                EData.Left(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
            end
        end
        for i=1:nright
            data_idx=depressed_right_pos(i)+ptsrange(1):depressed_right_pos(i)+ptsrange(2);
            for j=1:nCh
                EData.Right(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
            end
        end        
    case 'pressed'
        %data.EEG_Left{1}=zeros(nCh,ptsperCh,nleft);
        %data.EEG_Right{1}=zeros(nCh,ptsperCh,nright);
        EData.Left =zeros(ptsperCh,nCh,nleft);
        EData.Right=zeros(ptsperCh,nCh,nright);
        EData.Relax=[];
        if npress==0
            for i=1:nleft
                data_idx=pressed_left_pos(i)+ptsrange(1):pressed_left_pos(i)+ptsrange(2);
                for j=1:3
                    left_pressed_trial_idx=find(left_pressed_pos(j,:)==pressed_left_pos(i));
                    if ~isempty(left_pressed_trial_idx)
                        if left_outlier(left_pressed_trial_idx)
                            data_left_outlier(i)=true;
                            disp(['Skip: Outlier in ' num2str(i) ' of left pressed']);                        
                        end
                    end
                end
                for j=1:nCh                    
                    %data.EEG_Left{1}(j,:,i)=raweeg.EEG(IndChs(j),data_idx);
                    EData.Left(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
                end
            end
            for i=1:nright
                data_idx=pressed_right_pos(i)+ptsrange(1):pressed_right_pos(i)+ptsrange(2);
                for j=1:3
                    right_pressed_trial_idx=find(right_pressed_pos(j,:)==pressed_right_pos(i));
                    if ~isempty(right_pressed_trial_idx)
                        if right_outlier(right_pressed_trial_idx)
                            data_right_outlier(i)=true;
                            disp(['Skip: Outlier in ' num2str(i) ' of right pressed']);                        
                        end
                    end
                end
                for j=1:nCh
                    %data.EEG_Right{1}(j,:,i)=raweeg.EEG(IndChs(j),data_idx);
                    EData.Right(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
                end
            end
        else
            for i=1:nleft
                data_idx=left_pressed_pos(npress,i)+ptsrange(1):left_pressed_pos(npress,i)+ptsrange(2);
                if left_outlier(i)
                    data_left_outlier(i)=true;
                    disp(['Skip: Outlier in ' num2str(i) ' of left trial']);                        
                end                    
                for j=1:nCh
                    %data.EEG_Left{1}(j,:,i)=raweeg.EEG(IndChs(j),data_idx);
                    EData.Left(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
                end
            end
            for i=1:nright
                data_idx=right_pressed_pos(npress,i)+ptsrange(1):right_pressed_pos(npress,i)+ptsrange(2);
                if right_outlier(i)
                    data_right_outlier(i)=true;
                    disp(['Skip: Outlier in ' num2str(i) ' of right trial']);
                end
                for j=1:nCh
                    %data.EEG_Right{1}(j,:,i)=raweeg.EEG(IndChs(j),data_idx);
                    EData.Right(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
                end
            end            
        end
    case 'onesidepress' %for hemiplegia, one side is abled and "pressed" data to be extracted, the other is disabled and MI data to be extracted
        warning('In onesidepress mode, no outlier detection is applied.');
        if npress ~= 0 
            error('Currently the code doesn''t support multi-press in a single trial!\n');
        end
        TrialData = cell(1,2);
        event_pos = cell(1,2);
        iPosInTrial_MIEvent=1*raweeg.sampling_rate;
        if(~isempty(find(pressed_left_pos > 0, 1)) && ~isempty(find(pressed_right_pos > 0, 1))) % both classes have instances
            % find the class with less instances and use it as
            % the disabled side
            warning('Both sides have non-empty press records but you assume a hemiplegia? Anyway here we use the class with less instances as disabled side\n');
            if(length(find(pressed_left_pos>0)) < length(find(pressed_right_pos>0)) )
                iParalysisSide=1;iNonParalysisSide=2;
                event_pos{iParalysisSide}=indicate_left_pos(:)+iPosInTrial_MIEvent;
                event_pos{iNonParalysisSide}=pressed_right_pos;
            else
                iParalysisSide=2;iNonParalysisSide=1;
                event_pos{iParalysisSide}=indicate_right_pos(:)+iPosInTrial_MIEvent;
                event_pos{iNonParalysisSide}=pressed_left_pos;
            end
        else
            if(isempty(find(pressed_left_pos > 0, 1))) %left side paralysis
                iParalysisSide=1;iNonParalysisSide=2;
                event_pos{iParalysisSide}=indicate_left_pos(:)+iPosInTrial_MIEvent;
                event_pos{iNonParalysisSide}=pressed_right_pos;
            elseif(isempty(find(pressed_right_pos > 0, 1)))%right side paralysis
                iParalysisSide=2;iNonParalysisSide=1;
                event_pos{iParalysisSide}=indicate_right_pos(:)+iPosInTrial_MIEvent;
                event_pos{iNonParalysisSide}=pressed_left_pos;
            end
        end
        for iSide=1:2
            nTrial_iSide=length(event_pos{iSide});
            TrialData{iSide}=zeros(ptsperCh,nCh,nTrial_iSide);
            for i=1:nTrial_iSide
                data_idx=event_pos{iSide}(i)+ptsrange(1):event_pos{iSide}(i)+ptsrange(2);
                for j=1:nCh
                    TrialData{iSide}(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
                end
            end
        end
        if(iParalysisSide==1)
            EData.Left = TrialData{iParalysisSide};
            EData.Right = TrialData{iNonParalysisSide};
        else
            EData.Right = TrialData{iParalysisSide};
            EData.Left = TrialData{iNonParalysisSide};
        end
        EData.Relax=[];
    case 'indicate'
        %data.EEG_Left{1}=zeros(nCh,ptsperCh,nleft);
        %data.EEG_Right{1}=zeros(nCh,ptsperCh,nright);
        EData.Left=zeros(ptsperCh,nCh,nleft);
        EData.Right=zeros(ptsperCh,nCh,nright);
        EData.Relax=[];
        for i=1:nleft
            data_idx=indicate_left_pos(i)+ptsrange(1):indicate_left_pos(i)+ptsrange(2);
            if left_outlier(i)
                data_left_outlier(i)=true;
                disp(['Skip: Outlier in ' num2str(i) ' of left indicate']);                        
            end                    
            for j=1:nCh
                %data.EEG_Left{1}(j,:,i)=raweeg.EEG(IndChs(j),data_idx);
                EData.Left(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
                EData.Left(:,j,i)=[zeros(1,length(data_idx(data_idx<=0))) double(raweeg.EEG(IndChs(j),data_idx(data_idx>0)))*raweeg.resolution];
            end
        end
        for i=1:nright
            data_idx=indicate_right_pos(i)+ptsrange(1):indicate_right_pos(i)+ptsrange(2);
            if right_outlier(i)
                data_right_outlier(i)=true;
                disp(['Skip: Outlier in ' num2str(i) ' of right indicate']);
            end
            for j=1:nCh
                EData.Right(:,j,i)=[zeros(1,length(data_idx(data_idx<=0))) double(raweeg.EEG(IndChs(j),data_idx(data_idx>0)))*raweeg.resolution];
            end
        end
        if nrelax>0
            EData.Relax=zeros(ptsperCh,nCh,nrelax);
            for i=1:nrelax
                data_idx=indicate_relax_pos(i)+ptsrange(1):indicate_relax_pos(i)+ptsrange(2);
                % No outlier check for relaxed
                for j=1:nCh
                    %data.EEG_Right{1}(j,:,i)=raweeg.EEG(IndChs(j),data_idx);
                    EData.Relax(:,j,i)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
                end                
            end
        end
    case 'all'
        % Default extracts start of trial to 8s after start, 0.5s granularity
        % Each segment 4s long. Total 9 segments per trial.
        % data comprises left_data, left_event, right_data, right_event
        % data  = trials x segments x channels x time sample
        % event = 1 if key pressed, else, event=0
        Seg_length=4;
        Seg_gran=0.5;
        Seg_end_time=8;
        nSeg=(Seg_end_time-Seg_length)/Seg_gran+1;
        left_data=zeros(nleft,nSeg,nCh,Seg_length*raweeg.sampling_rate);
        left_event=zeros(nleft,nSeg);
        for i=1:nleft
            if left_outlier(i)
                data_left_outlier(i)=true;
                disp(['Skip: Outlier in ' num2str(i) ' of left indicate']);                        
            end
            for j=1:nSeg
                Seg_start=left_start_pos(i)+(j-1)*Seg_gran*raweeg.sampling_rate;
                Seg_end  =Seg_start+Seg_length*raweeg.sampling_rate-1;
                data_idx =Seg_start:Seg_end;
                last_seg_idx=Seg_end-(Seg_gran*raweeg.sampling_rate)+1:Seg_end;
                for k=1:3                   
                    if ~isempty(find(last_seg_idx==left_pressed_pos(k,i),1))
                        left_event(i,j)=1;
                    end
                end
                for k=1:nCh
                    left_data(i,j,k,:)=double(raweeg.EEG(IndChs(k),data_idx))*raweeg.resolution;
                end
            end
        end
        right_data=zeros(nright,nSeg,nCh,Seg_length*raweeg.sampling_rate);        
        right_event=zeros(nright,nSeg);
        for i=1:nright
            if right_outlier(i)
                data_right_outlier(i)=true;
                disp(['Skip: Outlier in ' num2str(i) ' of right indicate']);                        
            end
            for j=1:nSeg
                Seg_start=right_start_pos(i)+(j-1)*Seg_gran*raweeg.sampling_rate;
                Seg_end  =Seg_start+Seg_length*raweeg.sampling_rate-1;
                data_idx =Seg_start:Seg_end;
                for k=1:3                   
                    if ~isempty(find(last_seg_idx==right_pressed_pos(k,i),1))
                        right_event(i,j)=1;
                    end
                end
                for k=1:nCh
                    right_data(i,j,k,:)=double(raweeg.EEG(IndChs(k),data_idx))*raweeg.resolution;
                end
            end
        end
    case 'specific'
        nstim=length(specstim);
        extractpos = cell(nstim,1);
        nsamples=0;

        max_eeg_len = size(raweeg.EEG, 2); % ccwang: used to remove stim with not enough data
        
        for k=1:nstim
            extractpos{k} = raweeg.stimpos((raweeg.stimcode==specstim(k)));
            
            % ccwang: remove stims that is not long enough
            e_pos = extractpos{k} + ptsrange(2);
            rm_pts = find(e_pos > max_eeg_len);
            if ~isempty(rm_pts)
                fprintf('ccwang: stimcode %d removed because of length:', specstim(k));
                fprintf('%d ', rm_pts);
                fprintf('\n');
                extractpos{k}(rm_pts) = [];
            end
            clear e_pos;
            
            nsamples = nsamples+length(extractpos{k});
        end
        
        %data.x=zeros(ptsperCh,nCh,nsamples);
        data.x=zeros(nCh,ptsperCh,nsamples); %changed by ZY; will still get nSamples x nChs x nTrials in the end
        sidx=0;
        for k=1:nstim %for 3 stim codes
            for i=1:length(extractpos{k}) %for the nTrials in each stim code
                data_idx=extractpos{k}(i)+ptsrange(1):extractpos{k}(i)+ptsrange(2);
                sidx=sidx+1;
                % by default the data is nChs x nSamples x nTrials: added by ZY
                data.x(:,:,sidx)=double(raweeg.EEG(IndChs,data_idx)) * raweeg.resolution;                
%                 for j=1:nCh
%                     data.x(:,j,sidx)=double(raweeg.EEG(IndChs(j),data_idx))*raweeg.resolution;
%                 end
                data.y(sidx)=k-1;
            end
        end
        data.x = permute(data.x,[2 1 3]); %to convert to nSamples x nChs x nTrials; ZY
        data.c=strChList;        
        data.s=raweeg.sampling_rate;
end

% Remove outliers
switch xmethod
    case {'depressed','pressed','indicate'}
        %data.EEG_Left{1}(:,:,data_left_outlier)=[];
        %data.EEG_Right{1}(:,:,data_right_outlier)=[];
        EData.Left(:,:,data_left_outlier)=[];
        EData.Right(:,:,data_right_outlier)=[];
        if isempty(EData.Relax)
            data.x=cat(3,EData.Left,EData.Right);
            data.y=cat(2,zeros(1,size(EData.Left,3)),ones(1,size(EData.Right,3)));
        else
            data.x=cat(3,EData.Left,EData.Relax,EData.Right);
            data.y=cat(2,zeros(1,size(EData.Left,3)),0.5*ones(1,size(EData.Relax,3)),ones(1,size(EData.Right,3)));
        end
        data.c=strChList;        
        data.s=raweeg.sampling_rate;
    case{'onesidepress'}
        if isempty(EData.Relax)
            data.x=cat(3,EData.Left,EData.Right);
            data.y=cat(2,zeros(1,size(EData.Left,3)),ones(1,size(EData.Right,3)));
        else
            data.x=cat(3,EData.Left,EData.Relax,EData.Right);
            data.y=cat(2,zeros(1,size(EData.Left,3)),0.5*ones(1,size(EData.Relax,3)),ones(1,size(EData.Right,3)));
        end
        data.c=strChList;        
        data.s=raweeg.sampling_rate;
    case 'all'
        left_start_pos(data_left_outlier)=[];
        left_end_pos(data_left_outlier)=[];
        left_pressed_pos(:,data_left_outlier)=[];
        left_indicate_pos(data_left_outlier)=[];
        left_data(data_left_outlier,:,:,:)=[];
        left_event(data_left_outlier,:)=[];
        right_data(data_right_outlier,:,:,:)=[];
        right_event(data_right_outlier,:)=[];
        right_start_pos(data_right_outlier)=[];
        right_end_pos(data_right_outlier)=[];
        right_pressed_pos(:,data_right_outlier)=[];
        right_indicate_pos(data_right_outlier)=[];       
        data.EEG_Left=left_data;
        data.Left_Event=left_event;
        data.Left_Pos=[left_start_pos' left_indicate_pos' left_pressed_pos' left_end_pos'];      
        data.EEG_Right=right_data;
        data.Right_Event=right_event;
        data.Right_Pos=[right_start_pos' right_indicate_pos' right_pressed_pos' right_end_pos'];      
end

data.strSelChs = strChList; %raweeg.sel_chan_list;
data.sampling_rate = raweeg.sampling_rate;
% data.SampleTimes = ptsrange(1):ptsrange(2) / raweeg.sampling_rate;
