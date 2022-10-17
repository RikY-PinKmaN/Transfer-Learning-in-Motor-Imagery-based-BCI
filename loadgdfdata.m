function [data]=loadgdfdata(filename,ccode)
% LOADGDFDATA Loads data from BCI Competition 2a, to maintain consistency
% with existing data structure in I2R
%   by Chin Zheng Yang (zychin@i2r.a-star.edu.sg)
%   by Ang Kai Keng (kkang@i2r.a-star.edu.sg)
%
%   Syntax:
%     data=loadgdfdata(filename,ccode)
%   where
%     filename: GDF file with full path name
%     ccode:    769, 770, 771, 772
%     data:     Raw EEG data format
%
%   Example:
%     [Raw_BCI]=loadgdfdata('D:\KK\Data\BCI Competition 4\A01T.gdf',[769 770 771 772]);
%
%   check that you have installed BIOSIG first

bdebug = 0;
if bdebug == 1
    datadir = 'C:\Users\\goswa\Desktop\University\Peoject\gdf files';
    filename = 'C:\Users\goswa\Desktop\University\Peoject\gdf files\A01T.gdf';
   % C:\Users\Sidath\Documents\sourcecode\Mahnaz\BCICIV_2a_mat_NaNInterpolate
    ccode = [768 769 770 771 772];
end

rejcode = 1023;
nccodes = length(ccode);
if exist('sload','file')
    [s,h] = sload(filename, 0, 'OVERFLOWDETECTION:OFF');
else
    msgbox('please run biosig_installer in BIOSIG for sload to read in GDF files');
    error('please run biosig_installer in BIOSIG for sload to read in GDF files');
end
nch = size(s,2);
eventcodes = h.EVENT.TYP;
nevents = length(h.EVENT.TYP);
eventpos = h.EVENT.POS;
eventdur = h.EVENT.DUR;

%next get the stimpos of all trials (regardless valid or not)
idx = false(nevents,1);
codenum = zeros(1,nccodes);

for i=1:length(ccode)
    idx = idx |eventcodes == ccode(i);
end
%then get the stimcodes of all trials
stimcodes = eventcodes(idx);
stimpos = eventpos(idx);
stimdur = eventdur(idx);
BadTrialMarks = h.ArtifactSelection;
%now we only use the stimcodes of valid trials; removed to include all
%stimcodes = stimcodes(BadTrialMarks == 0);
%stimpos = stimpos(BadTrialMarks == 0);
for i=1:length(ccode)
    codenum(1,i)=sum(eventcodes == ccode(i));
end

data.numchannel=nch; %consistency with BCI Competition 3 Set 4a
data.num_channel=nch; %consistency with TEC data
data.sampling_rate=h.SampleRate;
data.resolution=1;%what is the resolution?? assume to be 1 (please confirm)
% BCI cnt data is uV x ch format. Raw EEG is ch x uV format, so require transpose
data.EEG=s';
data.stimpos=stimpos;
data.stimcode=stimcodes;
data.stimdur=stimdur;
data.codelist=ccode;
data.codenum=codenum;
%finally chan_list and sel_chan
switch nch
    case 25
        data.chan_list =         {'Fz';...
            'FC3';'FC1';'FCz';'FC2';'FC4';...
            'C5';'C3';'C1';'Cz';'C2';'C4';'C6';...
            'CP3';'CP1';'CPz';'CP2';'CP4';...
            'P1';'Pz';'P2';...
            'POz';'EOGL';'EOGC';'EOGR'};
        data.sel_chan_list =     {'Fz';...
            'FC3';'FC1';'FCz';'FC2';'FC4';...
            'C5';'C3';'C1';'Cz';'C2';'C4';'C6';...
            'CP3';'CP1';'CPz';'CP2';'CP4';...
            'P1';'Pz';'P2';...
            'POz';'EOGL';'EOGC';'EOGR'};
    case 6
        data.chan_list = {'C3';'Cz';'C4';'EOGL';'EOGC';'EOGR'};
        data.sel_chan_list = {'C3';'Cz';'C4';'EOGL';'EOGC';'EOGR'};
end
data.BadTrialMarks = BadTrialMarks;
data.bErrStimCode = false(1,length(stimcodes)); %due to no RemoveErrTrials


return