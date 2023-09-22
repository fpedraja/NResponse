%% read new files 3 (mimimc 1 2 and controls)
addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\6Fish_new_exp_2\data2\');
files2=dir([myKsDir, '\*NR_IDX_control*']);
files3=dir([myKsDir, '\*NR_IDX_test_1*']);
files4=dir([myKsDir, '\*NR_IDX_test_2*']);
c=0.55;

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


%% read feod for phth and plot

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\6Fish_new_exp\data2\');
files2=dir([myKsDir, '\*TIME_IDX_control*']);
files3=dir([myKsDir, '\*TIME_IDX_test_1*']);
files4=dir([myKsDir, '\*TIME_IDX_test_2*']);
c=0.65;

ANG=[]; NRTOT=[];
for i=1:size(files2,1)
   load([myKsDir,'\',files2(i).name]) 
   ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   
end
figure; 
% a=1;
% for i=1:180
%    try
%        plot(AUX(:,i),a,'.k'); hold on; a=a+1; 
%    end
% end

for t=1:32
    MEd=nanmedian(NRTOT(:,:,t),2);
    MEdN=1-(MEd./nanmean(MEd(1:20)));
     subplot(4,8,t); plot(-1:0.02:2,MEdN); ylim([-1 2]); hold on;
end
 


ANG=[]; NRTOT=[];
for i=1:size(files3,1)
   load([myKsDir,'\',files3(i).name]) 
   ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   
end

for t=1:32
    MEd=nanmedian(NRTOT(:,:,t),2);
    MEdN=1-(MEd./nanmean(MEd(1:20)));
     subplot(4,8,t); plot(-1:0.02:2,MEdN); ylim([-1 2]); hold on;
end    



ANG=[]; NRTOT=[];
for i=1:size(files4,1)
   load([myKsDir,'\',files4(i).name]) 
   ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   
end

for t=1:32
    MEd=nanmedian(NRTOT(:,:,t),2);
    MEdN=1-(MEd./nanmean(MEd(1:20)));
     subplot(4,8,t); plot(-1:0.02:2,MEdN); ylim([-1 2])
end   

%% fish change in pos head and tail all data day by day
addpath('D:\KIT3'); MEH1=[]; MEH2=[]; MEHC=[];
%%
%clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\6Fish_new_exp\');
files4=dir([myKsDir, '\video*.avi']);
load([myKsDir(1:end-9),'\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_EOD_data.mat'])

AUX2=unique(Obj_idx_1(:,2)); 
[Head_1,Tail_1]=head_tail_pos(myKsDir,files4,AUX2,Obj_idx_1);

AUX2=unique(Obj_idx_2(:,2)); 
[Head_2,Tail_2]=head_tail_pos(myKsDir,files4,AUX2,Obj_idx_2);

AUX2=unique(Obj_idx_control(:,2)); 
[Head_c,Tail_c]=head_tail_pos(myKsDir,files4,AUX2,Obj_idx_control);

mi=min([size(Head_1,1) size(Head_2,1) size(Head_c,1)]);
figure; boxplot([Head_1(1:mi,1) Head_2(1:mi,1) Head_c(1:mi,1)],'PlotStyle','compact')
meHead_1=median(Head_1); meHead_2=median(Head_2); meHead_c=median(Head_c); 
MEH1=[MEH1;meHead_1]; MEH2=[MEH2;meHead_2]; MEHC=[MEHC;meHead_c];

%%

figure; subplot(1,2,1); plot(MEH1,'-r'); hold on; plot(MEH2,'-r'); hold on; plot(MEHC,'-b');
subplot(1,2,2); boxplot([MEH1, MEH2, MEHC],'PlotStyle','compact')

figure; line([1 2 3],[MEH1, MEH2, MEHC])

AUX=zscore([MEH1, MEH2, MEHC],0,2); figure; boxplot(AUX,'PlotStyle','compact')

figure; line([1 2 3],[AUX])

