% %% Old analysis
% addpath('D:\KIT3');
% clearvars; %close all;
% myKsDir = uigetdir('Z:\locker\Fede\7Fish_new_exp_2\');
% files2=dir([myKsDir, '\EODdata2*']);
% files3=dir([myKsDir, '\obj_num*']);
% files4=dir([myKsDir, '\video*.avi']);
% load([myKsDir(1:end-9),'\Fish_5_',myKsDir(end-7:end),'_EOD_data.mat'])
% 
% 
% % for obj on and mimics on
% tic
% Ang=nan(30,16); NRe=nan(30,16); NRtotal=nan(30,16);
% a=1;
% for i=22:1:53
%     AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
%     [AUX1, ~]=find(Obj_idx(:,1)==i);
%     for j=1:size(AUX1,1)
%         [AUX2]=find(Obj_idx(:,2)==Obj_idx(AUX1(j),2));
%         M = readtable([myKsDir,'\CUT_' files4(Obj_idx(AUX2(1),2)).name(1:end-4) ,'DLC_resnet50_1Fish_2mimics_objMar18shuffle1_500000.csv']); M(1:2,:)=[];
%         
%         AUX4=Obj_idx(AUX2,1);
%         [AUX3]=find(AUX4==i); FrameIDX=((601*AUX3)-300);%place of video frame corresponding to the obj idx
%         fishH=table2array(M(FrameIDX-150:FrameIDX+150,2:3)); fishT=table2array(M(FrameIDX-150:FrameIDX+150,8:9));
%         head=[]; tail=[];
%         for t=1:301
%             head(t,1:2)=[str2num(fishH{t,1}) str2num(fishH{t,2})]; tail(t,1:2)=[str2num(fishT{t,1}) str2num(fishT{t,2})]; %figure; plot(head(1,1),head(1,2),'or',tail(1,1),tail(1,2),'.k'); %hold on; plot(tail())
%         end  
%         Head=median(head); Tail=median(tail);
%         Ang(j,a) = atan((Tail(1,2)-Head(1,2))/(Tail(1,1)-Head(1,1)));   
%         NRtotal(j,a)=-median(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+max([median(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=0 & Freqtotalpost(AUX1(j),:)<=1)), median(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=1 & Freqtotalpost(AUX1(j),:)<=2)) ]);
%        %AUX5=-mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=0 & Freqtotalpost(AUX1(j),:)<=1));
%        AUX5=-median(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+median(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=1 & Freqtotalpost(AUX1(j),:)<=2)); 
%        
%        if NRtotal(j,a)>0.25 || AUX5>0.25
%             NRe(j,a)=1;
%         else
%             NRe(j,a)=0;
%         end   
%     end
%     a=a+1;
% end
% 
% save(['Z:\locker\Fede\5Fish_new_exp\data2','\Fish_5_',myKsDir(end-7:end),'_NR_IDX_test.mat'],'Ang', 'NRe', 'NRtotal');
% 
% % for control data
% 
% Ang=nan(30,16); NRe=nan(30,16); NRtotal=nan(30,16);
% a=1;
% for i=22:2:52
%     AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
%     [AUX1, ~]=find(Obj_idx_control(:,1)==i);   
%     for j=1:size(AUX1,1)
%         [AUX2]=find(Obj_idx_control(:,2)==Obj_idx_control(AUX1(j),2));
%         M = readtable([myKsDir,'\CUT_' files4(Obj_idx_control(AUX2(1),2)).name(1:end-4) ,'DLC_resnet50_1Fish_2mimics_objMar18shuffle1_500000.csv']); M(1:2,:)=[];
%         
%         AUX4=Obj_idx_control(AUX2,1);
%         [AUX3]=find(AUX4==i); FrameIDX=((601*AUX3)-300);%place of video frame corresponding to the obj idx
%         fishH=table2array(M(FrameIDX-150:FrameIDX+150,2:3)); fishT=table2array(M(FrameIDX-150:FrameIDX+150,8:9));
%         head=[]; tail=[];
%         for t=1:301
%             head(t,1:2)=[str2num(fishH{t,1}) str2num(fishH{t,2})]; tail(t,1:2)=[str2num(fishT{t,1}) str2num(fishT{t,2})]; %figure; plot(head(1,1),head(1,2),'or',tail(1,1),tail(1,2),'.k'); %hold on; plot(tail())
%         end  
%         Head=median(head); Tail=median(tail);
%         Ang(j,a) = atan((Tail(1,2)-Head(1,2))/(Tail(1,1)-Head(1,1)));   
%         NRtotal(j,a)=-mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)<0 & Freqtotalpost__control(AUX1(j),:)>=-1))+mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)>=0 & Freqtotalpost__control(AUX1(j),:)<=1)); 
%         AUX5=-mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)<0 & Freqtotalpost__control(AUX1(j),:)>=-1))+mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)>=1 & Freqtotalpost__control(AUX1(j),:)<=2)); 
%        
%         if NRtotal(j,a)>0.25 || AUX5>0.25
%             NRe(j,a)=1;
%         else
%             NRe(j,a)=0;
%         end   
%     end
%     a=a+1;
% end
% 
% save(['Z:\locker\Fede\2Fish_mimic_obj\data2','\Fish_2_',myKsDir(end-7:end),'_NR_IDX_control.mat'],'Ang', 'NRe', 'NRtotal');
% toc
%% NEW IMPORTANT
addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\8Fish_new_exp_2\');
% FI=dir([myKsDir]);
% dirFlags = [FI.isdir];
% subDirs = FI(dirFlags);

files4=dir([myKsDir, '\video*.avi']);
load([myKsDir(1:end-9),'\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_EOD_data.mat'])


% for obj on and mimics 1 on
tic
Ang=nan(30,32); NRe=nan(30,32); NRtotal=nan(30,32);
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
Ang=nan(30,32); NRe=nan(30,32); NRtotal_1=nan(30,32);
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

Ang=nan(30,32); NRe=nan(30,32); NRtotal=nan(30,32);
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
%% OLD obj position check files from NR IDX angle
% addpath('D:\KIT3');
% clearvars; %close all;
% myKsDir = uigetdir('Z:\locker\Fede\2Fish_mimic_obj\');
% files2=dir([myKsDir, '\*NR_IDX_control*']);
% files3=dir([myKsDir, '\*NR_IDX_test_1*']);
% files4=dir([myKsDir, '\*NR_IDX_test_2*']);
% %c=0.52;
% 
% c=0.65; %0.65
% 
% POS=[1	1	26	1	26
% 2	1	23	1	23
% 3	1	20	1	20
% 4	1	17	1	17
% 5	1	14	1	14
% 6	1	11	1	11
% 7	1	8	1	8
% 8	1	5	1	5
% 9	5	1	5	1
% 10	8	1	8	1
% 11	11	1	11	1
% 12	14	1	14	1
% 13	17	1	17	1
% 14	20	1	20	1
% 15	23	1	23	1
% 16	26	1	26	1
% ];

%% check files from NR IDX angle
% addpath('D:\KIT3');
% clearvars; %close all;
% myKsDir = uigetdir('Z:\locker\Fede\2Fish_mimic_obj\');
% files2=dir([myKsDir, '\*NR_IDX_control*']);
% files3=dir([myKsDir, '\*NR_IDX_test*']);
% %c=0.52;
% 
% c=0.65; %0.65
% 
% POS=[1	16	1	1	16
% 2	11	1	1	11
% 3	6	1	1	6
% 4	1	1	1	1
% 5	1	6	6	1
% 6	1	11	11	1
% 7	1	16	16	1
% 8	4	4	4	4
% 9	4	9	9	4
% 10	4	14	14	4
% 11	9	4	4	9
% 12	9	9	9	9
% 13	14	4	4	14
% 14	6	11	11	6
% 15	14	14	14	14
% 16	16	16	16	16
% ];

%% control NR performance
ANG=[]; NRE=[]; NRTOT=[];
for i=1:size(files2,1)
   load([myKsDir,'\',files2(i).name]) 
   ANG=[ANG;Ang]; NRE=[ NRE; NRe]; NRTOT=[NRTOT; NRtotal];   
end

MA=nan(26,26,2); 
PA=nan(26,26,2);
for j=1:16 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRE(:,j); AUX(isnan(AUX))=[];
    AUX2=ANG(:,j); AUX2(isnan(AUX2))=[];
    
    AUX4=AUX(AUX2>0-c & AUX2<0+c); % x as x    
    AUX5=AUX(AUX2>1.5708-c | AUX2<-1.5708+c); % y as x
    
    MA(POS(j,2),POS(j,3),1)=sum(AUX4)/size(AUX4,1); 
    MA(POS(j,4),POS(j,5),2)=sum(AUX5)/size(AUX5,1);
    %MA(POS(j,2),POS(j,3),2)=sum(AUX5)/size(AUX5,1);
    
    PA(POS(j,2),POS(j,3),1)=size(AUX4,1); 
    PA(POS(j,4),POS(j,5),2)=size(AUX5,1);
    %PA(POS(j,2),POS(j,3),2)=size(AUX5,1); 
    
    
end

%% test NR performance

ANG=[]; NRE=[]; NRTOT=[];
for i=1:size(files3,1)
   load([myKsDir,'\',files3(i).name]) 
   ANG=[ANG;Ang]; NRE=[ NRE; NRe]; NRTOT=[NRTOT; NRtotal];   
end

MA_test=nan(16,16,2); 
PA_test=nan(16,16,2);

for j=1:16 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRE(:,j); AUX(isnan(AUX))=[];
    AUX2=ANG(:,j); AUX2(isnan(AUX2))=[];
    
    AUX4=AUX(AUX2>0-c & AUX2<0+c); % x as x    
    AUX5=AUX(AUX2>1.5708-c | AUX2<-1.5708+c); % y as x
    
    MA_test(POS(j,2),POS(j,3),1)=sum(AUX4)/size(AUX4,1);
    MA_test(POS(j,4),POS(j,5),2)=sum(AUX5)/size(AUX5,1);
    %MA_test(POS(j,2),POS(j,3),2)=sum(AUX5)/size(AUX5,1);
    
    PA_test(POS(j,2),POS(j,3),1)=size(AUX4,1); 
    PA_test(POS(j,4),POS(j,5),2)=size(AUX5,1);
    %PA_test(POS(j,2),POS(j,3),2)=size(AUX5,1);

end
%% plot from NR perfomance

SDS_1=inpaint_nans(nanmean(MA,3));
SDS2_1=inpaint_nans(nanmean(MA_test,3));

SDS_1=inpaint_nans(MA(:,:,2));
%SDS2_1=inpaint_nans(MA_test(:,:,2));

w     = 2;                              % Size of the sliding window (same number of cols and rows in this case)
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


figure(100); 
ax1=subplot(1,2,1); contourf(SDS,100,'edgecolor','none'); axis equal; colormap(ax1,brewermap([],'Blues')); %caxis([0 max(max(SDS))]);
ax2=subplot(1,2,2); contourf(SDS2,100,'edgecolor','none'); axis equal; colormap(ax2,brewermap([],'Reds')); %caxis([min() max(max(SDS2))]);

%% control zscore
ANG=[]; NRE=[]; NRTOT=[];
for i=1:size(files2,1)
   load([myKsDir,'\',files2(i).name]) 
   ANG=[ANG;Ang]; NRE=[ NRE; NRe]; NRTOT=[NRTOT; NRtotal];   
end

MA=nan(16,16,2); 

for j=1:16 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRTOT(:,j); AUX(isnan(AUX))=[];
    AUX2=ANG(:,j); AUX2(isnan(AUX2))=[];
    
    AUX4=AUX(AUX2>0-c & AUX2<0+c); % x as x    
    AUX5=AUX(AUX2>1.5708-c | AUX2<-1.5708+c); % y as x
    
    MA(POS(j,4),POS(j,5),2)=nanmean(AUX5);
    MA(POS(j,2),POS(j,3),1)=nanmean(AUX4);    
end

%% test zscore

ANG=[]; NRE=[]; NRTOT=[];
for i=1:size(files2,1)
   load([myKsDir,'\',files3(i).name]) 
   ANG=[ANG;Ang]; NRE=[ NRE; NRe]; NRTOT=[NRTOT; NRtotal];   
end

MA_test=nan(16,16,2); 

for j=1:16 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRTOT(:,j); AUX(isnan(AUX))=[];
    AUX2=ANG(:,j); AUX2(isnan(AUX2))=[];
    
    AUX4=AUX(AUX2>0-c & AUX2<0+c); % x as x    
    AUX5=AUX(AUX2>1.5708-c | AUX2<-1.5708+c); % y as x
    
    MA_test(POS(j,2),POS(j,3),1)=nanmean(AUX4);
    MA_test(POS(j,4),POS(j,5),2)=nanmean(AUX5);
end


%% plots from NR zscore

%SDS_1=inpaint_nans(nanmean(MA,3));
%SDS2_1=inpaint_nans(nanmean(MA_test,3));

SDS_1=inpaint_nans(MA(:,:,2));
SDS2_1=inpaint_nans(MA_test(:,:,2));

w = 2;                              % Size of the sliding window (same number of cols and rows in this case)
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


figure; 
ax1=subplot(1,2,1); contourf(SDS,100,'edgecolor','none'); axis equal; colormap(ax1,brewermap([],'Blues')); %caxis([0 max(max(SDS))]);
ax2=subplot(1,2,2); contourf(SDS2,100,'edgecolor','none'); axis equal; colormap(ax2,brewermap([],'Reds')); %caxis([min() max(max(SDS2))]);
