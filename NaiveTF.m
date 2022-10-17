%% Naive Transfer Learning
close all %Close any open figures
clc       %Clear command window
clear all %Clear workspace
for j=1:9
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
    subjectsi(j)=[];
    subjects(j)=[];
    test_subjects(j)=[];
  
    % 768 = start of trail, 769 = left arm cue onset
    % 770 = right arm cue onset 771 = foot cue onset
    % 772 = tongue cue onset 1023 = rejected trial
    allstim = [768 769 770 771 772 1023];
    stim = [769 770 771 772];

    dir2='C:\Users\goswa\Desktop\University\Project\true_labels'; %File path (change for your computer)
    datadir='C:\Users\goswa\Desktop\University\Peoject\matlab_files';
    timeseg=[0.5 2.5]; %Time segment to use (0.5 - 2.5 seconds after the cue).
    xR=[];
    xL=[];
    %% Source and target data
    for subi=1:8
        clear newdata; %Clear values at the beginning of each loop
        clear producedfeatur;
        clear V;

        %Load in the source data
        [Raw_sub,EMap]=loadbci4eegimagery([subjects{subi} '.mat'],[769 770]);
        xsubi=extracteegbci4imagery(Raw_sub,'indicate','seconds',[0 3],'selchs',EMap(1,subjectsi{subi}));
        xR=cat(3,xR,xsubi.Right(:,:,1:end));% Right Hand
        xL=cat(3,xL,xsubi.Left(:,:,1:end));% Left Hand
        clear Raw_sub;
    end
    %Load in the target data
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

    r=50;
    xR1=[];
    xL1=[];

    [Raw_sub,EMap]=loadbci4eegimagery([subjects{j} '.mat'],[769 770]);
    xsubi1=extracteegbci4imagery(Raw_sub,'indicate','seconds',[0 3],'selchs',EMap(1,subjectsi{j}));
    xR1=cat(3,xR1,xsubi1.Right(:,:,1:r)); % Right Hand
    xL1=cat(3,xL1,xsubi1.Left(:,:,1:r)); % Left Hand

    % Combining source and target data for training
    XR=cat(3,xR,xR1);
    XL=cat(3,xL,xL1);
    xsubi.x=cat(3,XR,XL);
    xsubi.y=[ones(size(XR,3),1);ones(size(XL,3),1)+1];
    % Testing data
    txsubi.x=cat(3,xsubi1.Right(:,:,r+1:end),xsubi1.Left(:,:,r+1:end));
    txsubi.y=[ones(size(xsubi1.Right(:,:,r+1:end),3),1);ones(size(xsubi1.Left(:,:,r+1:end),3),1)+1];
    clear Raw_sub; %Clear raw data

    ntrial=size(xR1,3); %Find the number of trials

    %% Bandpass filter the train and test data between 8 and 35Hz
    lowfreq(1,1)=8; %Lowest frequency to allow 8Hz
    highfreq(1,1)=35; %Highest frequency to allow 35Hz
    f=1;

    %Design a minimal order eliptcal filter with the frequencies shown
    %above and sampling rate of
    [n1,w1]=ellipord([lowfreq(1,f)/125,highfreq(1,f)/125],[(lowfreq(1,f)-1)/125,(highfreq(1,f)+1)/125],1,40);
    [b1,a1]=ellip(n1,1,40, w1); %return filter transfer function coefficients of filter oder n1

    %Filter to the train data (trn = train)
    finaltrn1.x=filter(b1,a1,xsubi.x);
    finaltrn.x=finaltrn1.x(126:625,:,:); %Use samples 125 to 625
    finaltrn.y=xsubi.y;

    %Filter test data
    finaltest.lowfreq=lowfreq(1,f); %Using the same bandpass filter as above (assuming f=1)
    finaltest.highfreq=highfreq(1,f);

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
    LDA_model1=fldt(trndataLDA); %Train the model using the train data
    LDA_model1=fldc(LDA_model1,tstdataLDA); %Test the model using test data
    TrainResult(1,j)=LDA_model1.taccuracy; %Compute training accuracy
    TestResult(1,j)=LDA_model1.caccuracy; %Compute testing accuracy
end
disp(mean(TestResult))