addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\');
files2=dir([myKsDir, '\EODdata2*']);

a=1; b=1; xs=[]; ys=[]; Freq1post=[]; Freq2post=[]; Freq1=[]; samplerate1=30000;  FreqTotal_1=[]; FreqTotal_2=[]; Freqtotalpost_1=[]; Freqtotalpost_2=[]; FreqTotal_control=[]; Freqtotalpost__control=[];
INTERVM2=[]; EODM2=[]; INTERVM1=[]; EODM1=[]; EODC=[];
%%
tic
for i=1:250%size(files2,1)
    fileID = fopen([myKsDir,'\',files2(i).name]);
    A = fread(fileID,[8,Inf],'int16'); fclose(fileID);
    [~,sample1_M1]=findpeaks(((A(1,:))) ,'MINPEAKHEIGHT',10000,'MINPEAKDISTANCE',100); % Mimic 1
    [~,sample1_M2]=findpeaks(((A(2,:))) ,'MINPEAKHEIGHT',10000,'MINPEAKDISTANCE',100); % Mimic 2
    if abs(min(zscore(A(6,:))))<=abs(max(zscore(A(6,:))))
        [~,sample1_FM1M2]=findpeaks((A(6,:)-min(A(6,:)))/ range(A(6,:)) ,'MINPEAKHEIGHT',0.55,'MINPEAKDISTANCE',100); % recording electrodes; fish + Mimic 1 + Mimic 2
    else
        [~,sample1_FM1M2]=findpeaks((-A(6,:)-min(-A(6,:)))/ range(A(6,:)) ,'MINPEAKHEIGHT',0.55,'MINPEAKDISTANCE',100);
    end
    clear A
    disp([i toc])
    %% move to no mimic or mimic part
    if (size(sample1_M1,2)<2000 || isempty(sample1_M1)==1) && (size(sample1_M2,2)<2000 || isempty(sample1_M2)==1)  % no mimic / control
        
        for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
        end
        
        for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
        end
        
        
        EODrate1=(diff(sample1_FM1M2)/samplerate1); %EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); %EODr1=zscore(EODr1);
        EODC=[EODC;EODrate1'];
        
    elseif (isempty(sample1_M2)==1 || size(sample1_M2,2)<500) && isempty(sample1_M1)==0     % mimic present at mimic 1
        for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
        end
        
        for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
        end
        EODrate1=(diff(sample1_FM1M2)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); %EODr1=zscore(EODr1);
        EODrateM1=(diff(sample1_M1)/samplerate1);
        
        a=1; Interv=[];
        for i=1:size(sample1_FM1M2,2)
            [idx2,idx]=min(abs(sample1_M1-sample1_FM1M2(i)));
            
            if sample1_M1(idx)-sample1_FM1M2(i)<0
                idx=idx+1;
            end
            try
                Interv(a)=(sample1_M1(idx)-sample1_FM1M2(i))/samplerate1;
                Interv(a+1)=(sample1_M1(idx-1)-sample1_FM1M2(i))/samplerate1;
                a=a+2;
            end
        end
        
        INTERVM1=[INTERVM1;Interv'];
        EODM1=[EODM1;EODrate1'];
        
    elseif (isempty(sample1_M1)==1 || size(sample1_M1,2)<500) && isempty(sample1_M2)==0    % mimic present at mimic 2
        for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
        end
        
        for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
        end
        
        EODrate1=(diff(sample1_FM1M2)/samplerate1); %EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); %EODr1=zscore(EODr1);
        EODrateM2=(diff(sample1_M2)/samplerate1);
        
        a=1; Interv=[];
        for i=1:size(sample1_FM1M2,2)
            [idx2,idx]=min(abs(sample1_M2-sample1_FM1M2(i)));
            
            if sample1_M2(idx)-sample1_FM1M2(i)<0
                idx=idx+1;
            end
            try
                Interv(a)=(sample1_M2(idx)-sample1_FM1M2(i))/samplerate1;
                Interv(a+1)=(sample1_M2(idx-1)-sample1_FM1M2(i))/samplerate1;
                a=a+2;
            end
        end
        
        INTERVM2=[INTERVM2;Interv'];
        EODM2=[EODM2;EODrate1'];
    end
end

EODC(EODC>0.5)=[]; EODM1(EODM1>0.5)=[]; EODM2(EODM2>0.5)=[];

figure; 
subplot(2,3,1); histogram(EODC,500) 
subplot(2,3,2); histogram(EODM1,500) 
subplot(2,3,3); histogram(EODM2,500) 
subplot(2,3,5); histogram(INTERVM1,500,'Normalization','pdf'); xlim([-0.05 0.05]); 
subplot(2,3,6); histogram(INTERVM2,500,'Normalization','pdf'); xlim([-0.05 0.05]); 

MC=1/median(EODC)
SEMC = std(EODC)/sqrt(length(EODC))
M1=1/median(EODM1)
SEM1 = std(EODM1)/sqrt(length(EODM1))
M2=1/median(EODM2)
SEM2 = std(EODM2)/sqrt(length(EODM2))

figure; 
subplot(2,3,1); histogram(1./EODC,500) 




