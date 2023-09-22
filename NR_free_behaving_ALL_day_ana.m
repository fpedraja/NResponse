addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\1Fish_mimic_obj\');
files2=dir([myKsDir, '\EODdata2*']);
%files3=dir([myKsDir, '\obj_num*']);
files4=dir([myKsDir, '\video*.avi']);
%%
a=0; b=1; xs=[]; ys=[]; samplerate1=30000; Time=[]; EOD=[]; 
figure; MEeod=[]; OP=[]; Meod=[];  M1total=[]; M2total=[];

for i=1:size(files2,1)
    tic
    fileID = fopen([myKsDir,'\',files2(i).name]);
    A = fread(fileID,[8,Inf],'int16'); fclose(fileID);
    [~,sample1_M1]=findpeaks(((A(1,:))) ,'MINPEAKHEIGHT',10000,'MINPEAKDISTANCE',100); % Mimic 1
    [~,sample1_M2]=findpeaks(((A(2,:))) ,'MINPEAKHEIGHT',10000,'MINPEAKDISTANCE',100); % Mimic 2
    [~,sample1_FM1M2]=findpeaks(((A(6,:))) ,'MINPEAKHEIGHT',11200,'MINPEAKDISTANCE',100);
    
    M1=[]; M2=[];
    
    if size(sample1_M1,2)<2000 || isempty(sample1_M1)==1% no mimic / control
        
        for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
        end
        
        for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
        end
        
        
        time=a+sample1_FM1M2/samplerate1;        
        a=a+size(A,2)/samplerate1; clear A
        EODrate1=(diff(sample1_FM1M2)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); %EODr1=zscore(EODr1);
        toc
       % plot(time,EODr1,'.b'); hold on; pause(.05);
        MEeod=[MEeod;nanmedian(EODr1)];
        Meod=[Meod;nanmean(EODr1)];
        OP=[OP;0];
    else % mimic present / real exp
        for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
        end
        
        for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
            sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
        end
        
        time=a+sample1_FM1M2/samplerate1;        
        a=a+size(A,2)/samplerate1; clear A
        EODrate1=(diff(sample1_FM1M2)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); %EODr1=zscore(EODr1);
        toc
        %plot(time,EODr1,'.r'); hold on; pause(.05);
        MEeod=[MEeod;nanmedian(EODr1)];
        Meod=[Meod;nanmean(EODr1)];
        OP=[OP;0];
        
        a=1;
        for k=1:size(sample1_FM1M2,2) %calculate IDI with fish and mimic
            try
%                 [idx2,idx]=min(abs(sample1_M1-sample1_FM1M2(k)));
%                 if sample1_M1(idx)>=sample1_FM1M2(k)
%                     M1(a,1)=(sample1_M1(idx-1)-sample1_FM1M2(k))/samplerate1;
%                     M1(a,2)=(sample1_M1(idx)-sample1_FM1M2(k))/samplerate1;
%                 else
%                     M1(a,1)=(sample1_M1(idx)-sample1_FM1M2(k))/samplerate1;
%                     M1(a,2)=(sample1_M1(idx+1)-sample1_FM1M2(k))/samplerate1;
%                 end
                
                [idx2,idx]=min(abs(sample1_M2-sample1_FM1M2(k)));
                if sample1_M2(idx)>=sample1_FM1M2(k)
                    M2(a,1)=(sample1_M2(idx-1)-sample1_FM1M2(k))/samplerate1;
                    M2(a,2)=(sample1_M2(idx)-sample1_FM1M2(k))/samplerate1;
                else
                    M2(a,1)=(sample1_M2(idx)-sample1_FM1M2(k))/samplerate1;
                    M2(a,2)=(sample1_M2(idx+1)-sample1_FM1M2(k))/samplerate1;
                end
                a=a+1;
            end
        end
        
    end
   
    Time=[Time,time];
    EOD=[EOD,EODr1];
    M1total=[M1total;M1];
    M2total=[M2total;M2];
      
end

figure; subplot(1,2,1); hist(M1total,500);
subplot(1,2,2); hist(M2total,500);

 











