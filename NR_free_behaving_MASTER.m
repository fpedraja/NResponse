% FIRST PART
addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\9Fish_new_exp_2\');
files2=dir([myKsDir, '\EODdata2*']);
files3=dir([myKsDir, '\obj_num*']);
files4=dir([myKsDir, '\video*.avi']);

a=1; b=1; xs=[]; ys=[]; Freq1post=[]; Freq2post=[]; Freq1=[]; samplerate1=30000;  FreqTotal_1=[]; FreqTotal_2=[]; Freqtotalpost_1=[]; Freqtotalpost_2=[]; FreqTotal_control=[]; Freqtotalpost__control=[];
Obj_idx_control=[]; Obj_idx_1=[]; Obj_idx_2=[];
%%
for i=1:size(files3,1)
    try %#ok<TRYNC> % check that the obj index has actual values (night vs day cycle)
        M = csvread([myKsDir,'\',files3(i).name]);
        if isempty(M)==1
        else
            fileID = fopen([myKsDir,'\',files2(i).name]);
            A = fread(fileID,[8,Inf],'int16'); fclose(fileID);
            [~,sample1_M1]=findpeaks(((A(1,:))) ,'MINPEAKHEIGHT',10000,'MINPEAKDISTANCE',100); % Mimic 1
            [~,sample1_M2]=findpeaks(((A(2,:))) ,'MINPEAKHEIGHT',10000,'MINPEAKDISTANCE',100); % Mimic 2
            [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)),'MINPEAKHEIGHT',0.8,'MINPEAKDISTANCE',60000); % find event, object-on
            
            
            if size(sample1_N,2)<size(M,1) %check if events in csv file corresponds to events detected
                [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)) ,'MINPEAKHEIGHT',0.7,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            elseif size(sample1_N,2)>size(M,1) %check if events in csv file corresponds to events detected
                [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)) ,'MINPEAKHEIGHT',0.9,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            end
            AUX23=0.5;
            while size(sample1_N,2) ~= size(M,1)
                if size(sample1_N,2)<size(M,1) %check if events in csv file corresponds to events detected
                    [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)) ,'MINPEAKHEIGHT',AUX23,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
                elseif size(sample1_N,2)>size(M,1) %check if events in csv file corresponds to events detected
                    [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)) ,'MINPEAKHEIGHT',AUX23,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
                end
                AUX23=AUX23+0.07;
                if AUX23 >=1
                    [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)) ,'MINPEAKHEIGHT',0.8,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
                    AUXN=sample1_N; sample1_N=[];
                    sample1_N=AUXN(2:end);
                    if size(sample1_N,2) ~= size(M,1)
                        AUXN=sample1_N; sample1_N=[];
                        sample1_N=AUXN(3:end);
                    end
                end
                if AUX23 >=1.05
                    break
                end
                
            end
            
            disp([size(sample1_N,2)-size(M,1) i])
            
            
            if  size(sample1_N,2)==size(M,1)
                if abs(min(zscore(A(6,:))))<=abs(max(zscore(A(6,:))))
                    [~,sample1_FM1M2]=findpeaks((A(6,:)-min(A(6,:)))/ range(A(6,:)) ,'MINPEAKHEIGHT',0.55,'MINPEAKDISTANCE',100); % recording electrodes; fish + Mimic 1 + Mimic 2
                else
                    [~,sample1_FM1M2]=findpeaks((-A(6,:)-min(-A(6,:)))/ range(A(6,:)) ,'MINPEAKHEIGHT',0.55,'MINPEAKDISTANCE',100);
                end
                %% move to no mimic or mimic part
                if (size(sample1_M1,2)<2000 || isempty(sample1_M1)==1) && (size(sample1_M2,2)<2000 || isempty(sample1_M2)==1)  % no mimic / control
                    
                    for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
                    end
                    
                    for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
                    end
                    
                    clear A
                    EODrate1=(diff(sample1_FM1M2)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); EODr1=zscore(EODr1);
                    
                    a=1; Freq1=[]; Freq1post=[]; AUX3=[]; AUX4=[];
                    vid = VideoReader([myKsDir,'\',files4(i).name]);
                    outputVideo = VideoWriter([myKsDir,'\CUT_',files4(i).name]);
                    outputVideo.FrameRate = vid.FrameRate;
                    open(outputVideo);
                    
                    for j=1:size(sample1_N,2)
                        [idx2,idx]=min(abs(sample1_FM1M2-sample1_N(j)));
                        try
                            Freq1post(a,:)=(sample1_FM1M2(idx-120:idx+220)-sample1_N(j))./samplerate1;
                            Freq1(a,:)=EODr1(idx-120:idx+220);
                            %[xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-60000/samplerate1 60000/samplerate1],0.9999999999);
                            %plot(xs,ys(a,:),'-k'); hold on;
                            for t=(round(sample1_N(j)/samplerate1*vid.FrameRate))-300:(round(sample1_N(j)/samplerate1*vid.FrameRate))+300
                                img = read(vid, t);  writeVideo(outputVideo,img);
                            end
                            AUX3(a,1)=M(j);
                            AUX4(a,1)=i;
                            a=a+1;
                        end
                    end
                    close(outputVideo);
                    FreqTotal_control=[FreqTotal_control;Freq1]; Freqtotalpost__control=[Freqtotalpost__control; Freq1post];
                    Obj_idx_control=[Obj_idx_control;[AUX3 AUX4]];
                    
                elseif (isempty(sample1_M2)==1 || size(sample1_M2,2)<500) && isempty(sample1_M1)==0     % mimic present at mimic 1
                    for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
                    end
                    
                    for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
                    end
                    
                    clear A
                    EODrate1=(diff(sample1_FM1M2)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); EODr1=zscore(EODr1);
                    
                    a=1; Freq1=[]; Freq1post_1=[]; AUX3=[]; AUX4=[];
                    vid = VideoReader([myKsDir,'\',files4(i).name]);
                    outputVideo = VideoWriter([myKsDir,'\CUT_',files4(i).name]);
                    outputVideo.FrameRate = vid.FrameRate;
                    open(outputVideo);
                    
                    for j=1:size(sample1_N,2)
                        [idx2,idx]=min(abs(sample1_FM1M2-sample1_N(j)));
                        try
                            Freq1post_1(a,:)=(sample1_FM1M2(idx-120:idx+220)-sample1_N(j))./samplerate1;
                            Freq1(a,:)=EODr1(idx-120:idx+220);
                            %[xs(a,:), ys(a,:)]=FitVal_EI(Freq1post_1(a,:), Freq1(a,:), [-60000/samplerate1 60000/samplerate1],0.9999999999);
                            %plot(xs,ys(a,:),'-k'); hold on;
                            for t=(round(sample1_N(j)/samplerate1*vid.FrameRate))-300:(round(sample1_N(j)/samplerate1*vid.FrameRate))+300
                                img = read(vid, t);  writeVideo(outputVideo,img);
                            end
                            AUX3(a,1)=M(j);
                            AUX4(a,1)=i;
                            a=a+1;
                        end
                    end
                    close(outputVideo);
                    FreqTotal_1=[FreqTotal_1;Freq1]; Freqtotalpost_1=[Freqtotalpost_1;Freq1post_1];
                    Obj_idx_1=[Obj_idx_1;[AUX3 AUX4]];
                    
                    
                elseif (isempty(sample1_M1)==1 || size(sample1_M1,2)<500) && isempty(sample1_M2)==0    % mimic present at mimic 2
                    for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
                    end
                    
                    for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
                    end
                    
                    clear A
                    EODrate1=(diff(sample1_FM1M2)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); EODr1=zscore(EODr1);
                    
                    a=1; Freq1=[]; Freq1post_2=[]; AUX3=[]; AUX4=[];
                    vid = VideoReader([myKsDir,'\',files4(i).name]);
                    outputVideo = VideoWriter([myKsDir,'\CUT_',files4(i).name]);
                    outputVideo.FrameRate = vid.FrameRate;
                    open(outputVideo);
                    
                    for j=1:size(sample1_N,2)
                        [idx2,idx]=min(abs(sample1_FM1M2-sample1_N(j)));
                        try
                            Freq1post_2(a,:)=(sample1_FM1M2(idx-120:idx+220)-sample1_N(j))./samplerate1;
                            Freq1(a,:)=EODr1(idx-120:idx+220);
                            % [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post_2(a,:), Freq1(a,:), [-60000/samplerate1 60000/samplerate1],0.9999999999);
                            %plot(xs,ys(a,:),'-k'); hold on;
                            for t=(round(sample1_N(j)/samplerate1*vid.FrameRate))-300:(round(sample1_N(j)/samplerate1*vid.FrameRate))+300
                                img = read(vid, t);  writeVideo(outputVideo,img);
                            end
                            AUX3(a,1)=M(j);
                            AUX4(a,1)=i;
                            a=a+1;
                        end
                    end
                    close(outputVideo);
                    FreqTotal_2=[FreqTotal_2;Freq1]; Freqtotalpost_2=[Freqtotalpost_2;Freq1post_2];
                    Obj_idx_2=[Obj_idx_2;[AUX3 AUX4]];
                    
                end               
            end
        end
    end
end

save([myKsDir(1:end-9),'\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_EOD_data.mat'],'FreqTotal_1','FreqTotal_2', 'FreqTotal_control', 'Freqtotalpost_1','Freqtotalpost_2', 'Freqtotalpost__control', 'Obj_idx_1','Obj_idx_2', 'Obj_idx_control')



% figure; subplot(2,2,1); plot(xs(1,:),ys,'-k');
% subplot(2,2,2);[hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');
% subplot(2,2,3);

figure; subplot(1,3,1)
for t=1:size(Freqtotalpost_1,1)
    plot(Freqtotalpost_1(t,:),t,'.k'); xlim([-1 2]); hold on;
end

subplot(1,3,2)
for t=1:size(Freqtotalpost_2,1)
    plot(Freqtotalpost_2(t,:),t,'.k'); xlim([-1 2]); hold on;
end

subplot(1,3,3)
for t=1:size(Freqtotalpost__control,1)
    plot(Freqtotalpost__control(t,:),t,'.k'); xlim([-1 2]); hold on;
end
% subplot(2,2,4);
% for i=1:size(Freq2post,1)
%     plot(Freq2post(i,:),i,'.k'); xlim([-1 2]); hold on;
% end
%% SECOND PART

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\9Fish_new_exp_2\');
% FI=dir([myKsDir]);
% dirFlags = [FI.isdir];
% subDirs = FI(dirFlags);

files4=dir([myKsDir, '\video*.avi']);
load([myKsDir(1:end-9),'\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_EOD_data.mat'])


% for obj on and mimics 1 on
tic
Ang=nan(50,32); NRe=nan(50,32); NRtotal=nan(50,32);
a=1;
for i=22:1:53
    AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
    [AUX1, ~]=find(Obj_idx_1(:,1)==i);
    for j=1:size(AUX1,1)
        [AUX2]=find(Obj_idx_1(:,2)==Obj_idx_1(AUX1(j),2));
        M = readtable([myKsDir,'\CUT_' files4(Obj_idx_1(AUX2(1),2)).name(1:end-4) ,'DLC_resnet50_NEW_expe_3May11shuffle1_300000.csv']); M(1:2,:)=[];
        
        AUX4=Obj_idx_1(AUX2,1);
        [AUX3]=find(AUX4==i); FrameIDX=((601*AUX3)-300);%place of video frame corresponding to the obj idx
        fishH=table2array(M(FrameIDX-150:FrameIDX+150,2:3)); fishT=table2array(M(FrameIDX-150:FrameIDX+150,8:9));
        head=[]; tail=[];
        for t=1:301
            head(t,1:2)=[str2num(fishH{t,1}) str2num(fishH{t,2})]; tail(t,1:2)=[str2num(fishT{t,1}) str2num(fishT{t,2})]; %figure; plot(head(1,1),head(1,2),'or',tail(1,1),tail(1,2),'.k'); %hold on; plot(tail())
        end  
        Head=median(head); Tail=median(tail);
        Ang(j,a) = atan2((Tail(1,2)-Head(1,2)),(Tail(1,1)-Head(1,1)));   
        NRtotal(j,a)=-mean(FreqTotal_1(AUX1(j),Freqtotalpost_1(AUX1(j),:)<0 & Freqtotalpost_1(AUX1(j),:)>=-1))+max([mean(FreqTotal_1(AUX1(j),Freqtotalpost_1(AUX1(j),:)>=0 & Freqtotalpost_1(AUX1(j),:)<=1)) mean(FreqTotal_1(AUX1(j),Freqtotalpost_1(AUX1(j),:)>=1 & Freqtotalpost_1(AUX1(j),:)<=2))]);
       %AUX5=-mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=0 & Freqtotalpost(AUX1(j),:)<=1));
       AUX5=-mean(FreqTotal_1(AUX1(j),Freqtotalpost_1(AUX1(j),:)<0 & Freqtotalpost_1(AUX1(j),:)>=-1))+mean(FreqTotal_1(AUX1(j),Freqtotalpost_1(AUX1(j),:)>=1 & Freqtotalpost_1(AUX1(j),:)<=2)); 
       
       if NRtotal(j,a)>0.20 || AUX5>0.20
            NRe(j,a)=1;
        else
            NRe(j,a)=0;
        end   
    end
    a=a+1;
end

save([myKsDir(1:end-9),'\data2','\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_NR_IDX_test_1.mat'],'Ang', 'NRe', 'NRtotal');

% for object on and mimic 2 on
Ang=nan(50,32); NRe=nan(50,32); NRtotal_1=nan(50,32);
a=1;
for i=22:1:53
    AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
    [AUX1, ~]=find(Obj_idx_2(:,1)==i);
    for j=1:size(AUX1,1)
        [AUX2]=find(Obj_idx_2(:,2)==Obj_idx_2(AUX1(j),2));
        M = readtable([myKsDir,'\CUT_' files4(Obj_idx_2(AUX2(1),2)).name(1:end-4) ,'DLC_resnet50_NEW_expe_3May11shuffle1_300000.csv']); M(1:2,:)=[];
        
        AUX4=Obj_idx_2(AUX2,1);
        [AUX3]=find(AUX4==i); FrameIDX=((601*AUX3)-300);%place of video frame corresponding to the obj idx
        fishH=table2array(M(FrameIDX-150:FrameIDX+150,2:3)); fishT=table2array(M(FrameIDX-150:FrameIDX+150,8:9));
        head=[]; tail=[];
        for t=1:301
            head(t,1:2)=[str2num(fishH{t,1}) str2num(fishH{t,2})]; tail(t,1:2)=[str2num(fishT{t,1}) str2num(fishT{t,2})]; %figure; plot(head(1,1),head(1,2),'or',tail(1,1),tail(1,2),'.k'); %hold on; plot(tail())
        end  
        Head=median(head); Tail=median(tail);
        Ang(j,a) = atan2((Tail(1,2)-Head(1,2)),(Tail(1,1)-Head(1,1)));   
        NRtotal(j,a)=-mean(FreqTotal_2(AUX1(j),Freqtotalpost_2(AUX1(j),:)<0 & Freqtotalpost_2(AUX1(j),:)>=-1))+max([mean(FreqTotal_2(AUX1(j),Freqtotalpost_2(AUX1(j),:)>=0 & Freqtotalpost_2(AUX1(j),:)<=1)) mean(FreqTotal_2(AUX1(j),Freqtotalpost_2(AUX1(j),:)>=1 & Freqtotalpost_2(AUX1(j),:)<=2))]);
       %AUX5=-mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=0 & Freqtotalpost(AUX1(j),:)<=1));
       AUX5=-mean(FreqTotal_2(AUX1(j),Freqtotalpost_2(AUX1(j),:)<0 & Freqtotalpost_2(AUX1(j),:)>=-1))+mean(FreqTotal_2(AUX1(j),Freqtotalpost_2(AUX1(j),:)>=1 & Freqtotalpost_2(AUX1(j),:)<=2)); 
       
       if NRtotal(j,a)>0.20 || AUX5>0.20
            NRe(j,a)=1;
        else
            NRe(j,a)=0;
        end   
    end
    a=a+1;
end

save([myKsDir(1:end-9),'\data2','\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_NR_IDX_test_2.mat'],'Ang', 'NRe', 'NRtotal');


% for control data

Ang=nan(50,32); NRe=nan(50,32); NRtotal=nan(50,32);
a=1;
for i=22:1:53
    AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
    [AUX1, ~]=find(Obj_idx_control(:,1)==i);   
    for j=1:size(AUX1,1)
        [AUX2]=find(Obj_idx_control(:,2)==Obj_idx_control(AUX1(j),2));
        M = readtable([myKsDir,'\CUT_' files4(Obj_idx_control(AUX2(1),2)).name(1:end-4) ,'DLC_resnet50_NEW_expe_3May11shuffle1_300000.csv']); M(1:2,:)=[];
        
        AUX4=Obj_idx_control(AUX2,1);
        [AUX3]=find(AUX4==i); FrameIDX=((601*AUX3)-300);%place of video frame corresponding to the obj idx
        fishH=table2array(M(FrameIDX-150:FrameIDX+150,2:3)); fishT=table2array(M(FrameIDX-150:FrameIDX+150,8:9));
        head=[]; tail=[];
        for t=1:301
            head(t,1:2)=[str2num(fishH{t,1}) str2num(fishH{t,2})]; tail(t,1:2)=[str2num(fishT{t,1}) str2num(fishT{t,2})]; %figure; plot(head(1,1),head(1,2),'or',tail(1,1),tail(1,2),'.k'); %hold on; plot(tail())
        end  
        Head=median(head); Tail=median(tail);
        Ang(j,a) = atan2((Tail(1,2)-Head(1,2)),(Tail(1,1)-Head(1,1)));   
        NRtotal(j,a)=-mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)<0 & Freqtotalpost__control(AUX1(j),:)>=-1))+max([mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)>=0 & Freqtotalpost__control(AUX1(j),:)<=1)) mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)>=1 & Freqtotalpost__control(AUX1(j),:)<=2))]); 
        AUX5=-mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)<0 & Freqtotalpost__control(AUX1(j),:)>=-1))+mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)>=1 & Freqtotalpost__control(AUX1(j),:)<=2)); 
       
        if NRtotal(j,a)>0.20 || AUX5>0.20
            NRe(j,a)=1;
        else
            NRe(j,a)=0;
        end   
    end
    a=a+1;
end

save([myKsDir(1:end-9),'\data2','\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_NR_IDX_control.mat'],'Ang', 'NRe', 'NRtotal');
toc

%% THIRD PART
%% read new files 3 (mimimc 1 2 and controls)
addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\9Fish_new_exp_2\data2\');
files2=dir([myKsDir, '\*NR_IDX_control*']);
files3=dir([myKsDir, '\*NR_IDX_test_1*']);
files4=dir([myKsDir, '\*NR_IDX_test_2*']);
c=0.55; % angle range for data exclusion 

POS=[1	26	1	26	1
2	25	2	25	2
3	23	1	23	1
4	22	2	22	2
5	20	1	20	1
6	19	2	19	2
7	17	1	17	1
8	16	2	16	2
9	14	1	14	1
10	13	2	13	2
11	11	1	11	1
12	10	2	10	2
13	8	1	8	1
14	7	2	7	2
15	5	1	5	1
16	4	2	4	2
17	1	4	1	4
18	2	5	2	5
19	1	7	1	7
20	2	8	2	8
21	1	10	1	10
22	2	11	2	11
23	1	13	1	13
24	2	14	2	14
25	1	16	1	16
26	2	17	2	17
27	1	19	1	19
28	2	20	2	20
29	1	22	1	22
30	2	23	2	23
31	1	25	1	25
32	2	26	2	26
];

figure;
%% mimc 2
ANG=[]; NRE=[]; NRTOT=[];
for i=1:size(files4,1)
   load([myKsDir,'\',files4(i).name]) 
   ANG=[ANG;Ang]; NRE=[ NRE; NRe]; NRTOT=[NRTOT; NRtotal];   
end

MA_test2=nan(26,26,2); 
PA_test2=nan(26,26,2);

for j=1:32 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRE(:,j);
    AUX2=ANG(:,j); AUX2(isnan(AUX))=[];  AUX(isnan(AUX))=[];
    
    AUX4=AUX(AUX2>0-c & AUX2<0+c); % x as x    
    AUX5=AUX(AUX2>1.5708-c | AUX2<-1.5708+c); % y as x
    
    MA_test2(POS(j,2),POS(j,3),1)=sum(AUX4)/size(AUX4,1); 
    MA_test2(POS(j,4),POS(j,5),2)=sum(AUX5)/size(AUX5,1);
    %MA(POS(j,2),POS(j,3),2)=sum(AUX5)/size(AUX5,1);
    
    PA_test2(POS(j,2),POS(j,3),1)=size(AUX4,1); 
    PA_test2(POS(j,4),POS(j,5),2)=size(AUX5,1);
    %PA(POS(j,2),POS(j,3),2)=size(AUX5,1); 
     
end

AUX2 = reshape(ANG,[1,size(ANG,1)*size(ANG,2)]); AUX2(isnan(AUX2))=[]; %AUX2(AUX2>1.5708-c | AUX2<-1.5708+c)=[];
addpath('C:\Users\fedu1\Google Drive\CircStat2012a')
subAx4 = subplot(1,3,3, polaraxes);
obj2 = CircHist(rad2deg(AUX2), 360,'parent', subAx4);
obj2.colorBar = 'k';  % change color of bars
obj2.avgAngH.LineStyle = '--'; % make average-angle line dashed
obj2.avgAngH.LineWidth = 1; % make average-angle line thinner
obj2.colorAvgAng = [.5 .5 .5]; % change average-angle line color
obj2.polarAxs.ThetaZeroLocation = 'bottom';
rl = rlim; % get current limits
delete(obj2.rH)
obj2.drawArrow(obj2.avgAng, obj2.r * range(rl), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'r')
obj2.drawScale; % update scale bar

%% mimic 1
ANG=[]; NRE=[]; NRTOT=[];
for i=1:size(files3,1)
   load([myKsDir,'\',files3(i).name]) 
   ANG=[ANG;Ang]; NRE=[ NRE; NRe]; NRTOT=[NRTOT; NRtotal];   
end

MA_test1=nan(26,26,2); 
PA_test1=nan(26,26,2);

for j=1:32 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRE(:,j); 
    AUX2=ANG(:,j); AUX2(isnan(AUX))=[]; AUX(isnan(AUX))=[];
    
    AUX4=AUX(AUX2>0-c & AUX2<0+c); % x as x    
    AUX5=AUX(AUX2>1.5708-c | AUX2<-1.5708+c); % y as x
    
    MA_test1(POS(j,2),POS(j,3),1)=sum(AUX4)/size(AUX4,1); 
    MA_test1(POS(j,4),POS(j,5),2)=sum(AUX5)/size(AUX5,1);
    %MA(POS(j,2),POS(j,3),2)=sum(AUX5)/size(AUX5,1);
    
    PA_test1(POS(j,2),POS(j,3),1)=size(AUX4,1); 
    PA_test1(POS(j,4),POS(j,5),2)=size(AUX5,1);
    %PA(POS(j,2),POS(j,3),2)=size(AUX5,1); 
     
end

AUX2 = reshape(ANG,[1,size(ANG,1)*size(ANG,2)]); AUX2(isnan(AUX2))=[]; %AUX2(AUX2>1.5708-c | AUX2<-1.5708+c)=[];
addpath('C:\Users\fedu1\Google Drive\CircStat2012a')
subAx4 = subplot(1,3,2, polaraxes);
obj2 = CircHist(rad2deg(AUX2), 360,'parent', subAx4);
obj2.colorBar = 'k';  % change color of bars
obj2.avgAngH.LineStyle = '--'; % make average-angle line dashed
obj2.avgAngH.LineWidth = 1; % make average-angle line thinner
obj2.colorAvgAng = [.5 .5 .5]; % change average-angle line color
obj2.polarAxs.ThetaZeroLocation = 'bottom';
rl = rlim; % get current limits
delete(obj2.rH)
obj2.drawArrow(obj2.avgAng, obj2.r * range(rl), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'r')
obj2.drawScale; % update scale bar
%% control
ANG=[]; NRE=[]; NRTOT=[];
for i=1:size(files2,1)
   load([myKsDir,'\',files2(i).name]) 
   ANG=[ANG;Ang]; NRE=[ NRE; NRe]; NRTOT=[NRTOT; NRtotal];   
end

MA_control=nan(26,26,2); 
PA_control=nan(26,26,2);

for j=1:32 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRE(:,j); 
    AUX2=ANG(:,j); AUX2(isnan(AUX))=[]; AUX(isnan(AUX))=[];
    
    AUX4=AUX(AUX2>0-c & AUX2<0+c); % x as x    
    AUX5=AUX(AUX2>1.5708-c | AUX2<-1.5708+c); % y as x
    
    MA_control(POS(j,2),POS(j,3),1)=sum(AUX4)/size(AUX4,1); 
    MA_control(POS(j,4),POS(j,5),2)=sum(AUX5)/size(AUX5,1);
    %MA(POS(j,2),POS(j,3),2)=sum(AUX5)/size(AUX5,1);
    
    PA_control(POS(j,2),POS(j,3),1)=size(AUX4,1); 
    PA_control(POS(j,4),POS(j,5),2)=size(AUX5,1);
    %PA(POS(j,2),POS(j,3),2)=size(AUX5,1); 
     
end

AUX2 = reshape(ANG,[1,size(ANG,1)*size(ANG,2)]); AUX2(isnan(AUX2))=[]; %AUX2(AUX2>1.5708-c | AUX2<-1.5708+c)=[];
addpath('C:\Users\fedu1\Google Drive\CircStat2012a')
subAx4 = subplot(1,3,1, polaraxes);
obj2 = CircHist(rad2deg(AUX2), 360,'parent', subAx4);
obj2.colorBar = 'k';  % change color of bars
obj2.avgAngH.LineStyle = '--'; % make average-angle line dashed
obj2.avgAngH.LineWidth = 1; % make average-angle line thinner
obj2.colorAvgAng = [.5 .5 .5]; % change average-angle line color
obj2.polarAxs.ThetaZeroLocation = 'bottom';
rl = rlim; % get current limits
delete(obj2.rH)
obj2.drawArrow(obj2.avgAng, obj2.r * range(rl), 'HeadWidth', 10, 'LineWidth', 2, 'Color', 'r')
obj2.drawScale; % update scale bar



%% interpolate
SDS_1=inpaint_nans(MA_test1(:,:,2));
SDS2_1=inpaint_nans(MA_test2(:,:,2));
SDS3_1=inpaint_nans(MA_control(:,:,2));

w     = 5;   % Size of the sliding window (same number of cols and rows in this case)
% Extrapolate values for current window
[Nr,Nc] = size(SDS_1);
Nextra  = 0.5*(w-1);
Ap      = interp2(1:Nc,1:Nr,SDS_1,-Nextra+1:Nc+Nextra,(-Nextra+1:Nr+Nextra).','makima');    % 2D extrapolation must use 'spline' or 'makima' interpolation
% Smooth data with sliding window
H  = ones(w)./w^2;                      % The 2D averaging filter
SDS  = filter2(H,Ap,'valid'); 

Ap      = interp2(1:Nc,1:Nr,SDS2_1,-Nextra+1:Nc+Nextra,(-Nextra+1:Nr+Nextra).','makima');    % 2D extrapolation must use 'spline' or 'makima' interpolation
SDS2  = filter2(H,Ap,'valid'); 
% The smooth resulting matrix

Ap      = interp2(1:Nc,1:Nr,SDS3_1,-Nextra+1:Nc+Nextra,(-Nextra+1:Nr+Nextra).','makima');    % 2D extrapolation must use 'spline' or 'makima' interpolation
SDS3  = filter2(H,Ap,'valid'); 

% figure; 
% ax1=subplot(3,1,1); contourf(xbins(1:end-1),ybins(1:end-1),sds,100,'edgecolor','none'); axis equal; colormap(ax1,brewermap([],'Blues')); %caxis([0 max(max(SDS))]);
% ax2=subplot(3,1,2); contourf(xbins(1:end-1),ybins(1:end-1),sds2,100,'edgecolor','none'); axis equal; colormap(ax2,brewermap([],'Reds')); %caxis([min() max(max(SDS2))]);
% ax3=subplot(3,1,3); contourf(xbins(1:end-1),ybins(1:end-1),sds3,100,'edgecolor','none'); axis equal; colormap(ax3,brewermap([],'Greys')); %caxis([0 max(max(SDS3))]);

% xbins=1:27;
% ybins=1:27;
% figure(100); 
% ax1=subplot(1,3,1); contourf(xbins(1:end-1),ybins(1:end-1),SDS,100,'edgecolor','none'); axis equal; colormap(ax1,brewermap([],'Reds')); %caxis([0 max(max(SDS))]);
% ax2=subplot(1,3,2); contourf(xbins(1:end-1),ybins(1:end-1),SDS2,100,'edgecolor','none'); axis equal; colormap(ax2,brewermap([],'Reds')); %caxis([min() max(max(SDS2))]);
% ax3=subplot(1,3,3); contourf(xbins(1:end-1),ybins(1:end-1),SDS3,100,'edgecolor','none'); axis equal; colormap(ax3,brewermap([],'Blues')); %caxis([0 max(max(SDS3))]);



MA_test1b=nan(26,26); MA_test2b=nan(26,26); MA_testcontrolb=nan(26,26); 
for t=1:size(POS,1)
MA_test1b(POS(t,2),POS(t,3))=SDS(POS(t,2),POS(t,3));
MA_test2b(POS(t,2),POS(t,3))=SDS2(POS(t,2),POS(t,3));
MA_testcontrolb(POS(t,2),POS(t,3))=SDS3(POS(t,2),POS(t,3));
end


%% FOURTH PART plots

%% plot for test 1 and 2
Irate=MA_test1b(:,:,1);
IrateNUM=PA_test1(:,:,2);

x=1:26;y=1:26;
colormap2 = brewermap(sum(~isnan(reshape(Irate,1,[]))),'Greys');alpha=0.6;alpha2=0.6;
s=sort(Irate(~isnan(reshape(Irate,1,[]))));circle=10;
numPoints=100;
theta=linspace(0,2*pi,numPoints); %100 evenly spaced points between 0 and 2pi
% imshow(cat(3,Hintergrund,Hintergrund,Hintergrund)); hold on
figure;
for i=1:size(x,2)
    for j=1:size(y,2)
        if ~isnan(Irate(i,j))
            col=colormap2((find(Irate(i,j)<=s,1)),:);
            fPos = get(gcf, 'Position');
            xl = xlim(); yl = ylim();
%             w = circle*(xl(2)-xl(1))/fPos(3);
%             h = circle*(yl(2)-yl(1))/fPos(4);

            w = circle*IrateNUM(i,j)/fPos(3);
            h = circle*IrateNUM(i,j)/fPos(4);
            
            mx = h*sin(theta); my = h*cos(theta);
            patch(x(i)+mx*2, y(j)+my*2,col, 'FaceColor', col, 'EdgeColor', 'none');
        else
            col=colormap2((find(Irate(i,j)<=s,1)),:);
            fPos = get(gcf, 'Position');
            xl = xlim(); yl = ylim();
            w = 2*(xl(2)-xl(1))/fPos(3);h = 2*(yl(2)-yl(1))/fPos(4);
            mx = w*sin(theta); my = h*cos(theta);
            patch(x(i)+mx, y(j)+my*0.8,[0 0 0], 'FaceColor', [0 0 0],  'EdgeColor', [0 0 0]);
        end
    end
end
hcb=colorbar; axis equal
set(get(colorbar,'ylabel'),'String', 'EOD rate')

%% plots for controls

Irate=MA_testcontrolb(:,:,1);
  IrateNUM=PA_control(:,:,2);

x=1:26;y=1:26;
colormap2 = brewermap(sum(~isnan(reshape(Irate,1,[]))),'Greys');alpha=0.6;alpha2=0.6;
s=sort(Irate(~isnan(reshape(Irate,1,[]))));circle=10;
numPoints=100;
theta=linspace(0,2*pi,numPoints); %100 evenly spaced points between 0 and 2pi
% imshow(cat(3,Hintergrund,Hintergrund,Hintergrund)); hold on
figure;
for i=1:size(x,2)
    for j=1:size(y,2)
        if ~isnan(Irate(i,j))
            col=colormap2((find(Irate(i,j)<=s,1)),:);
            fPos = get(gcf, 'Position');
            xl = xlim(); yl = ylim();
%             w = circle*(xl(2)-xl(1))/fPos(3);
%             h = circle*(yl(2)-yl(1))/fPos(4);

            w = circle*IrateNUM(i,j)/fPos(3);
            h = circle*IrateNUM(i,j)/fPos(4);
            
            mx = h*sin(theta); my = h*cos(theta);
            patch(x(i)+mx*2, y(j)+my*2,col, 'FaceColor', col, 'EdgeColor', 'none');
        else
            col=colormap2((find(Irate(i,j)<=s,1)),:);
            fPos = get(gcf, 'Position');
            xl = xlim(); yl = ylim();
            w = 2*(xl(2)-xl(1))/fPos(3);h = 2*(yl(2)-yl(1))/fPos(4);
            mx = w*sin(theta); my = h*cos(theta);
            patch(x(i)+mx, y(j)+my*0.8,[0 0 0], 'FaceColor', [0 0 0],  'EdgeColor', [0 0 0]);
        end
    end
end
hcb=colorbar; axis equal
set(get(colorbar,'ylabel'),'String', 'EOD rate')