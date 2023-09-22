%% to run raster plots and phth of eod f
addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\8Fish_new_exp\');
files4=dir([myKsDir, '\video*.avi']);
load([myKsDir(1:end-9),'\Fish_' myKsDir(16) '_',myKsDir(end-7:end),'_EOD_data.mat'])


unique(Obj_idx_1);

%
% for obj on and mimics 1 on
tic
Ang=nan(30,32);
Time=nan(181,30,32); EODtime=nan(341,30,32);
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
        Ang(j,a) = atan((Tail(1,2)-Head(1,2))/(Tail(1,1)-Head(1,1)));
        [xs, ys]=FitVal_EI(Freqtotalpost_1(AUX1(j),:),FreqTotal_1(AUX1(j),:),  [-1.3 2.3],0.99);
        Time(:,j,a)=ys;
        EODtime(:,j,a)=Freqtotalpost_1(AUX1(j),:);
    end
    a=a+1;
end

save([myKsDir(1:end-9),'\data2','\Fish_' myKsDir(16) '_',myKsDir(end-7:end),'_TIME_IDX_test_1.mat'],'Ang', 'Time', 'EODtime');

% for object on and mimic 2 on
Ang=nan(30,32);
Time=nan(181,30,32); EODtime=nan(341,30,32);
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
        Ang(j,a) = atan((Tail(1,2)-Head(1,2))/(Tail(1,1)-Head(1,1)));
        [xs, ys]=FitVal_EI(Freqtotalpost_2(AUX1(j),:),FreqTotal_2(AUX1(j),:),  [-1.3 2.3],0.99);
        Time(:,j,a)=ys;
        EODtime(:,j,a)=Freqtotalpost_2(AUX1(j),:);
        
    end
    a=a+1;
end

save([myKsDir(1:end-9),'\data2','\Fish_'  myKsDir(16) '_',myKsDir(end-7:end),'_TIME_IDX_test_2.mat'],'Ang', 'Time', 'EODtime');


% for control data

Ang=nan(30,32);
Time=nan(181,30,32); EODtime=nan(341,30,32);
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
        Ang(j,a) = atan((Tail(1,2)-Head(1,2))/(Tail(1,1)-Head(1,1)));
        [xs, ys]=FitVal_EI(Freqtotalpost__control(AUX1(j),:),FreqTotal_control(AUX1(j),:),  [-1.3 2.3],0.99);
        Time(:,j,a)=ys;
        EODtime(:,j,a)=Freqtotalpost__control(AUX1(j),:);
        
    end
    a=a+1;
end

save([myKsDir(1:end-9),'\data2','\Fish_'  myKsDir(16) '_',myKsDir(end-7:end),'_TIME_IDX_control.mat'],'Ang', 'Time', 'EODtime');
toc


%% read feod for phth and plot

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\6Fish_new_exp\data2\');
files2=dir([myKsDir, '\*TIME_IDX_control*']);
files3=dir([myKsDir, '\*TIME_IDX_test_1*']);
files4=dir([myKsDir, '\*TIME_IDX_test_2*']);
c=0.65;

%%
for j=15:19
    figure;
    ANG=[]; NRTOT=[]; EOD=[];
    for i=1:size(files2,1)
        load([myKsDir,'\',files2(i).name])
        ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   EOD=[EOD EODtime];
    end

    a=1; AUX=EOD(:,:,j);
    for i=1:180
        if isnan(AUX(:,i))==0
            subplot(2,3,1); plot(AUX(:,i),a,'.k'); hold on; a=a+1; xlim([-1 2])
        end
    end
    
    % for t=1:32
        MEd=nanmedian(NRTOT(:,:,j),2);
        MEdN=1-(MEd./nanmean(MEd(1:20)));
        Mad=nanstd(NRTOT(:,:,j),[],2);
         subplot(2,3,4); [hl, hp]=boundedline(-1.3:0.02:2.3, MEdN,Mad,'-k'); %plot(-1:0.02:2,MEdN); ylim([-1 2]); hold on;
    % end
    %
    
    ANG=[]; NRTOT=[]; EOD=[];
    for i=1:size(files3,1)
        load([myKsDir,'\',files3(i).name])
        ANG=[ANG;Ang]; NRTOT=[NRTOT Time];    EOD=[EOD EODtime];
    end
    a=1; AUX=EOD(:,:,j);
    for i=1:180
        if isnan(AUX(:,i))==0
            subplot(2,3,2); plot(AUX(:,i),a,'.k'); hold on; a=a+1; xlim([-1 2])
        end
    end
    
    % for t=1:32
        MEd=nanmedian(NRTOT(:,:,j),2);
        Mad=nanstd(NRTOT(:,:,j),[],2);
        MEdN=1-(MEd./nanmean(MEd(1:20)));
         subplot(2,3,5); [hl, hp]=boundedline(-1.3:0.02:2.3, MEdN,Mad,'-k'); %plot(-1:0.02:2,MEdN); ylim([-1 2]); hold on;
    % end
    
    ANG=[]; NRTOT=[]; EOD=[];
    for i=1:size(files4,1)
        load([myKsDir,'\',files4(i).name])
        ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   EOD=[EOD EODtime];
    end
    a=1; AUX=EOD(:,:,j);
    for i=1:180
        if isnan(AUX(:,i))==0
            subplot(2,3,3); plot(AUX(:,i),a,'.k'); hold on; a=a+1; xlim([-1 2])
        end
    end
    
    % for t=1:32
        MEd=nanmedian(NRTOT(:,:,j),2);
        Mad=nanstd(NRTOT(:,:,j),[],2);
        MEdN=1-(MEd./nanmean(MEd(1:20)));
         subplot(2,3,6); [hl, hp]=boundedline(-1.3:0.02:2.3, MEdN,Mad,'-k'); %plot(-1:0.02:2,MEdN); ylim([-1 2])
    % end
end