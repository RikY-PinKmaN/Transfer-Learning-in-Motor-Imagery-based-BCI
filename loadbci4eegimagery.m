function varargout=loadbci4eegimagery(filename,ccode)
load(filename);
% nch= number of channel
nch=size(raweeg.EEG,1);
npos=length(raweeg.trial_cue_pos);
%  mrk.y=mrk.y-min(mrk.y)+1;
% if isempty(find(mrk.y==2,1)) && ~isempty(find(mrk.y==3,1))
%     mrk.y(mrk.y==3)=2;
% end
varargout{1}.numchannels=nch;
varargout{1}.sampling_rate=raweeg.sampling_rate;
varargout{1}.resolution=raweeg.resolution;
varargout{1}.EEG=raweeg.EEG;
varargout{1}.stimpos=raweeg.trial_cue_pos;
varargout{1}.trial_start_pos=raweeg.trial_start_pos;
varargout{1}.trial_end_pos=raweeg.trial_end_pos;
varargout{1}.trial_class=raweeg.trial_class;
varargout{1}.trial_num_nan=raweeg.trial_num_nan;
varargout{1}.trial_nan_EOG=raweeg.trial_nan_EOG;
varargout{1}.trial_bad_mark=raweeg.trial_bad_mark;
varargout{1}.trial_num_blink=raweeg.trial_num_blink;
varargout{1}.eye_start_pos=raweeg.eye_start_pos;
varargout{1}.eye_end_pos=raweeg.eye_end_pos;
varargout{1}.eye_num_blink=raweeg.eye_num_blink;
varargout{1}.eye_num_nan=raweeg.eye_num_nan;
varargout{1}.eye_nan_EOG=raweeg.eye_nan_EOG;
varargout{1}.eye_cue_label=raweeg.eye_cue_label;
varargout{1}.nan_interp=raweeg.nan_interp;
varargout{1}.trial_cue_pos=raweeg.trial_cue_pos;
% Convert to stimcode
varargout{1}.stimcode=zeros(1,length(ccode));
varargout{1}.codelist=ccode;
varargout{1}.codenum=zeros(1,length(ccode));
for i=1:length(ccode)
    varargout{1}.stimcode(raweeg.trial_class==i)=ccode(i);
    varargout{1}.codenum(i)=length(find(raweeg.trial_class==i));
end
varargout{1}.chan_list=raweeg.chan_list';
varargout{1}.sel_chan_list=raweeg.chan_list';
if nargout>1
    varargout{2}=cell(3,nch);
    varargout{2}(1,:)=raweeg.chan_list';
    xpos=[0,-0.3629,-0.1821,0,0.1821,0.3629,-0.5769,-0.3846,-0.1923,0,0.1923,0.3846,0.5769,-0.3629,-0.1821,0,0.1821,0.3629,-0.1503,0,0.1503,0,0,0,0];
   varargout{2}(2,:)=num2cell(xpos);
    varargout{2}(3,:)=num2cell([0.3846,0.2126,0.1971,0.1923,0.1971,0.2126,0,0,0,0,0,0,0,-0.2126,-0.1971,-0.1923,-0.1971,-0.2126,-0.3928,-0.3846,-0.3928,-0.5769,0,0,0]);
end
