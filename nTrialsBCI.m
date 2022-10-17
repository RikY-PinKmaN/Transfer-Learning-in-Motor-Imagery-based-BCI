%% State of the art classifer with n_trials
clc       %Clear command window
clear all %Clear workspace

%Subject names
subjects={'A01T','A02T','A03T','A04T','A05T','A06T','A07T','A08T','A09T'};
test_subjects={'A01E','A02E','A03E','A04E','A05E','A06E','A07E','A08E','A09E'};
A01T=unique([1:22]);
A02T=unique([1:22]);
A03T=unique([1:22]);
A04T=unique([1:22]);
A05T=unique([1:22]);
A06T=unique([1:22]);
A07T=unique([1:22]);
A08T=unique([1:22]);
A09T=unique([1:22]);
subjectsi={A01T,A02T,A03T,A04T,A05T,A06T,A07T,A08T,A09T};

% 768 = start of trail, 769 = left arm cue onset
% 770 = right arm cue onset 771 = foot cue onset
% 772 = tongue cue onset 1023 = rejected trial
allstim = [768 769 770 771 772 1023];
stim = [769 770 771 772];

dir2='C:\Users\goswa\Desktop\University\Project\true_labels'; %File path (change for your computer)
datadir='C:\Users\goswa\Desktop\University\Peoject\matlab_files';
timeseg=[0.5 2.5]; %Time segment to use (0.5 - 2.5 seconds after the cue).

%% Data for each of the 9 subjetcs
for subi=1:9
    clear newdata; %Clear values at the beginning of each loop
    clear producedfeatur;
    clear V;
    r=5;
    %Load in the data
    [Raw_sub,EMap]=loadbci4eegimagery([subjects{subi} '.mat'],[769 770]); % RIght and Left Hand MI
    xsubi1=extracteegbci4imagery(Raw_sub,'indicate','seconds',[0 3],'selchs',EMap(1,subjectsi{subi}));
    % Training data
    xsubi.x=cat(3,xsubi1.Right(:,:,1:r),xsubi1.Left(:,:,1:r));
    xsubi.y=[ones(size(xsubi1.Right(:,:,1:r),3),1);ones(size(xsubi1.Left(:,:,1:r),3),1)+1];
    % Testing Data
    txsubi.x=cat(3,xsubi1.Right(:,:,r+1:end),xsubi1.Left(:,:,r+1:end));
    txsubi.y=[ones(size(xsubi1.Right(:,:,r+1:end),3),1);ones(size(xsubi1.Left(:,:,r+1:end),3),1)+1];

    clear Raw_sub; %Clear raw data
    ntrial=size(xsubi.x,3); %Find the number of trials
    f=1;
    %% Bandpass filter the train and test data between 8 and 35Hz
    lowfreq(1,1)=8; %Lowest frequency to allow 8Hz
    highfreq(1,1)=35; %Highest frequency to allow 35Hz

    %Design a minimal order eliptcal filter with the frequencies shown
    %above and sampling rate of
    [n1,w1]=ellipord([lowfreq(1,f)/125,highfreq(1,f)/125],[(lowfreq(1,f)-1)/125,(highfreq(1,f)+1)/125],1,40);
    [b1,a1]=ellip(n1,1,40, w1); %return filter transfer function coefficients of filter oder n1

    %Filter to the train data (trn = train)
    finaltrn1.x=filter(b1,a1,xsubi.x);
    finaltrn.x=finaltrn1.x(126:625,:,:); %Use samples 125 to 625
    finaltrn.y=xsubi.y;

    finaltest.lowfreq=lowfreq(1,f); %Using the same bandpass filter as above (assuming f=1)
    finaltest.highfreq=highfreq(1,f);

    %Filter test data
    finaltest1.x=filter(b1,a1,txsubi.x);
    finaltest.x=finaltest1.x(126:625,:,:); %Use samples 125 to 625
    finaltest.y=txsubi.y(:,:);

    %% Common spatial patterns and linear disciminant analysis machine learning
    numchannel=22;
    %Common spatial patterns (CSP) feature selection
    [finalfeaturtrn,selectedw1]=featcrossval(finaltrn,[1:22],numchannel); %training features
    finalfeaturtest=featcrostest(finaltest,[1:22],numchannel,selectedw1); %test features

    trndataLDA=[finalfeaturtrn.x' finalfeaturtrn.y-1]; %train data
    tstdataLDA=[finalfeaturtest.x' finalfeaturtest.y-1]; %Test data

    % Train and classify using fisher linear discriminnat analysis
    LDA_model1{subi}=fldt(trndataLDA); %Train the model using the train data
    LDA_model1{subi}=fldc(LDA_model1{subi},tstdataLDA); %Test the model using test data
    TrainResult(1,subi)=LDA_model1{subi}.taccuracy; %Compute training accuracy
    TestResult(1,subi)=LDA_model1{subi}.caccuracy; %Compute testing accuracy

end
% Average accuracy
disp(mean(TestResult))