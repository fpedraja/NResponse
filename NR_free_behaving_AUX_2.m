%% tail and head distribution
addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\1Fish_mimic_obj\');
files2=dir([myKsDir, '\EODdata2*']);
files3=dir([myKsDir, '\obj_num*']);
files4=dir([myKsDir, '\video*.avi']);
load(['Z:\locker\Fede\1Fish_mimic_obj','\Fish_5_',myKsDir(end-7:end),'_EOD_data.mat'])

%% for obj on and mimics on
tic
Ang=nan(30,16); NRe=nan(30,16); NRtotal=nan(30,16); TAIL=[]; HEAD=[];
a=1; b=1;
for i=22:2:52
    AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
    [AUX1, ~]=find(Obj_idx(:,1)==i);
    for j=1:size(AUX1,1)
        [AUX2]=find(Obj_idx(:,2)==Obj_idx(AUX1(j),2));
        M = readtable([myKsDir,'\CUT_' files4(Obj_idx(AUX2(1),2)).name(1:end-4) ,'DLC_resnet50_1Fish_2mimics_objMar18shuffle1_500000.csv']); M(1:2,:)=[];
        
        AUX4=Obj_idx(AUX2,1);
        [AUX3]=find(AUX4==i); FrameIDX=((601*AUX3)-300);%place of video frame corresponding to the obj idx
        fishH=table2array(M(FrameIDX-150:FrameIDX+150,2:3)); fishT=table2array(M(FrameIDX-150:FrameIDX+150,8:9));
        head=[]; tail=[];
        for t=1:301
            head(t,1:2)=[str2num(fishH{t,1}) str2num(fishH{t,2})]; tail(t,1:2)=[str2num(fishT{t,1}) str2num(fishT{t,2})]; %figure; plot(head(1,1),head(1,2),'or',tail(1,1),tail(1,2),'.k'); %hold on; plot(tail())
        end  
        Head=median(head); Tail=median(tail);
        Ang(j,a) = atan((Tail(1,2)-Head(1,2))/(Tail(1,1)-Head(1,1)));   
        NRtotal(j,a)=-mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=0 & Freqtotalpost(AUX1(j),:)<=1));
       %AUX5=-mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=0 & Freqtotalpost(AUX1(j),:)<=1));
       AUX5=-mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=1 & Freqtotalpost(AUX1(j),:)<=2)); 
       
       if NRtotal(j,a)>0.25 || AUX5>0.25
            NRe(j,a)=1;
        else
            NRe(j,a)=0;
       end
       
       TAIL(b,:)=Tail(1,:);
       HEAD(b,:)=Head(1,:);
        b=b+1;
    end
    a=a+1;
end

save(['Z:\locker\Fede\1Fish_mimic_obj\data2','\Fish_5_',myKsDir(end-7:end),'_HEAD_TAIL_IDX_test.mat'],'Ang', 'NRe', 'NRtotal','TAIL','HEAD');

% for control data

Ang=nan(30,16); NRe=nan(30,16); NRtotal=nan(30,16); TAIL=[]; HEAD=[];
a=1; b=1;
for i=22:2:52
    AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
    [AUX1, ~]=find(Obj_idx_control(:,1)==i);   
    for j=1:size(AUX1,1)
        [AUX2]=find(Obj_idx_control(:,2)==Obj_idx_control(AUX1(j),2));
        M = readtable([myKsDir,'\CUT_' files4(Obj_idx_control(AUX2(1),2)).name(1:end-4) ,'DLC_resnet50_1Fish_2mimics_objMar18shuffle1_500000.csv']); M(1:2,:)=[];
        
        AUX4=Obj_idx_control(AUX2,1);
        [AUX3]=find(AUX4==i); FrameIDX=((601*AUX3)-300);%place of video frame corresponding to the obj idx
        fishH=table2array(M(FrameIDX-150:FrameIDX+150,2:3)); fishT=table2array(M(FrameIDX-150:FrameIDX+150,8:9));
        head=[]; tail=[];
        for t=1:301
            head(t,1:2)=[str2num(fishH{t,1}) str2num(fishH{t,2})]; tail(t,1:2)=[str2num(fishT{t,1}) str2num(fishT{t,2})]; %figure; plot(head(1,1),head(1,2),'or',tail(1,1),tail(1,2),'.k'); %hold on; plot(tail())
        end  
        Head=median(head); Tail=median(tail);
        Ang(j,a) = atan((Tail(1,2)-Head(1,2))/(Tail(1,1)-Head(1,1)));   
        NRtotal(j,a)=-mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)<0 & Freqtotalpost__control(AUX1(j),:)>=-1))+mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)>=0 & Freqtotalpost__control(AUX1(j),:)<=1)); 
        AUX5=-mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)<0 & Freqtotalpost__control(AUX1(j),:)>=-1))+mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)>=1 & Freqtotalpost__control(AUX1(j),:)<=2)); 
       
        if NRtotal(j,a)>0.25 || AUX5>0.25
            NRe(j,a)=1;
        else
            NRe(j,a)=0;
        end 
        
       TAIL(b,:)=Tail(1,:);
       HEAD(b,:)=Head(1,:);
        b=b+1;
    end
    a=a+1;
end

save(['Z:\locker\Fede\1Fish_mimic_obj\data2','\Fish_5_',myKsDir(end-7:end),'_HEAD_TAIL_IDX_control.mat'],'Ang', 'NRe', 'NRtotal','TAIL','HEAD');
toc

%% read head and tail files

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\1Fish_mimic_obj\');
files2=dir([myKsDir, '\*_HEAD_TAIL_IDX_control*']);
files3=dir([myKsDir, '\*_HEAD_TAIL_IDX_test*']);
%% put all head and tail files together
Tail=[]; Head=[];
for i=1:size(files2,1)
   load([myKsDir,'\',files2(i).name]) 
   Head=[Head;HEAD]; Tail=[ Tail; TAIL];  
end
%% run gaussinan fit to head data
Data=Head;
X = 0:0.5:752;
Y = 0:0.5:832;

D = length(Data(1,:));
Mu = mean(Data);
Sigma = cov(Data);
P_Gaussian = zeros(length(X),length(Y)); AUX=zeros(length(X),length(Y));
for i=1:length(X)
   for j1=1:length(Y)
       x = [X(i),Y(j1)];
       P_Gaussian(i,j1) = 1/((2*pi)^(D/2)*sqrt(det(Sigma)))...
                    *exp(-1/2*(x-Mu)*Sigma^-1*(x-Mu)');
     AUX(i,j1)=sum(Data(:,1)>=X(i) & Data(:,1)<X(i)+1 & Data(:,2)>=Y(j1) & Data(:,2)<Y(j1)+1);
   end
end
figure; contourb(Y,X,P_Gaussian);
xlim([0 832]); ylim([0 758]); axis equal; 

figure; contourb(Y,X,AUX); 
xlim([0 832]); ylim([0 758]); axis equal; 


%% run gaussinan fit to tail data
Data=Tail;
X = 0:10:752;
Y = 0:10:832;

D = length(Data(1,:));
Mu = mean(Data);
Sigma = cov(Data);
P_Gaussian = zeros(length(X),length(Y)); AUX=zeros(length(X),length(Y));
for i=1:length(X)-1
   for j1=1:length(Y)-1
       x = [X(i),Y(j1)];
       P_Gaussian(i,j1) = 1/((2*pi)^(D/2)*sqrt(det(Sigma)))...
                    *exp(-1/2*(x-Mu)*Sigma^-1*(x-Mu)');
       AUX(i,j1)=sum(Data(:,1)>=X(i) & Data(:,1)<X(i+1) & Data(:,2)>=Y(j1) & Data(:,2)<Y(j1+1));
   end
end
figure; contourb(P_Gaussian,10); 
xlim([0 832]); ylim([0 758]); axis equal; 

figure; contourb(AUX); 
xlim([0 832]); ylim([0 758]); axis equal; 


figure; subplot(2,2,1); [N,c]=hist3(Head,'Nbins',[300 300],'CDataMode','auto','FaceColor','interp');
subplot(2,2,2); [N,c]=hist3(Tail,'Nbins',[300 300],'CDataMode','auto','FaceColor','interp')

%% test
Tail=[]; Head=[];
for i=1:size(files3,1)
   load([myKsDir,'\',files3(i).name]) 
   Head=[Head;HEAD]; Tail=[ Tail; TAIL];  
end

Data=Head;
X = 0:1:752;
Y = 0:1:832;

D = length(Data(1,:));
Mu = mean(Data);
Sigma = cov(Data);
P_Gaussian = zeros(length(X),length(Y)); AUX=zeros(length(X),length(Y));
for i=1:length(X)
   for j1=1:length(Y)
       x = [X(i),Y(j1)];
       P_Gaussian(i,j1) = 1/((2*pi)^(D/2)*sqrt(det(Sigma)))...
                    *exp(-1/2*(x-Mu)*Sigma^-1*(x-Mu)');
     AUX(i,j1)=sum(Data(:,1)>=X(i) & Data(:,1)<X(i)+1 & Data(:,2)>=Y(j1) & Data(:,2)<Y(j1)+1);
   end
end
figure; contourb(P_Gaussian,100);
xlim([0 832]); ylim([0 758]); axis equal; 

figure; contourb(AUX); 
xlim([0 832]); ylim([0 758]); axis equal; 


%%
Data=Tail(:,1:2);
X = 0:1:752;
Y = 0:1:832;

D = length(Data(1,:));
Mu = mean(Data);
Sigma = cov(Data);
P_Gaussian = zeros(length(X),length(Y)); AUX=zeros(length(X),length(Y));
for i=1:length(X)
   for j1=1:length(Y)
       x = [X(i),Y(j1)];
       P_Gaussian(i,j1) = 1/((2*pi)^(D/2)*sqrt(det(Sigma)))...
                    *exp(-1/2*(x-Mu)*Sigma^-1*(x-Mu)');
       AUX(i,j1)=sum(Data(:,1)>=X(i) & Data(:,1)<X(i)+1 & Data(:,2)>=Y(j1) & Data(:,2)<Y(j1)+1);
   end
end
figure; contourb(P_Gaussian,100); 
xlim([0 832]); ylim([0 758]); axis equal; 

figure; contourb(AUX); 
xlim([0 832]); ylim([0 758]); axis equal; 


subplot(2,2,3); hist3(Head,'Nbins',[30 30],'CDataMode','auto','FaceColor','interp')
subplot(2,2,4); hist3(Tail(:,1:2),'Nbins',[30 30],'CDataMode','auto','FaceColor','interp')
%%
figure;
for i=1:size(Tail,1)
    line([Head(i,1);Tail(i,1)],[Head(i,2);Tail(i,2)],'LineWidth',1.75,'Color',[0 0.4470 0.7410 0.1] )
    hold on;
    plot(Head(i,1),Head(i,2),'.k','MarkerSize',10,'Color',[0 0 0 0.1]); axis equal
end
    
%% head and tail run NEW IMPORTANT

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\7Fish_new_exp\');
files4=dir([myKsDir, '\video*.avi']);
load([myKsDir(1:end-9),'\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_EOD_data.mat'])


% for obj on and mimics 1 on
tic
a=1; b=1; TAIL=[]; HEAD=[]; FISHSP=[];
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
        FISHSP(b,1)=abs(sum(median(head(1:150,:))-median(head(152:end,:))))+abs(sum(median(tail(1:150,:))-median(tail(152:end,:))));
        TAIL(b,:)=[Tail(1,:) i];
        HEAD(b,:)=Head(1,:);
        b=b+1;
    end
    a=a+1;
end

save([myKsDir(1:end-9),'\data2','\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_HEAD_IDX_test_1.mat'],'HEAD', 'TAIL','FISHSP');

% for object on and mimic 2 on
a=1; b=1; TAIL=[]; HEAD=[];
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
        FISHSP(b,1)=abs(sum(median(head(1:150,:))-median(head(152:end,:))))+abs(sum(median(tail(1:150,:))-median(tail(152:end,:))));
        TAIL(b,:)=[Tail(1,:) i];
        HEAD(b,:)=Head(1,:);
        b=b+1;
        
    end
    a=a+1;
end

save([myKsDir(1:end-9),'\data2','\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_HEAD_IDX_test_2.mat'],'HEAD', 'TAIL','FISHSP');


% for control data

a=1; b=1; TAIL=[]; HEAD=[];
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
        FISHSP(b,1)=abs(sum(median(head(1:150,:))-median(head(152:end,:))))+abs(sum(median(tail(1:150,:))-median(tail(152:end,:))));
        TAIL(b,:)=[Tail(1,:) i];
        HEAD(b,:)=Head(1,:);
        b=b+1;
    end
    a=a+1;
end

save([myKsDir(1:end-9),'\data2','\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_HEAD_IDX_control.mat'],'HEAD', 'TAIL','FISHSP');
toc

%% NEW head ana

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\5Fish_new_exp\');
files2=dir([myKsDir, '\*_HEAD_IDX_control*']);
files3=dir([myKsDir, '\*_HEAD_IDX_test_1*']);
files4=dir([myKsDir, '\*_HEAD_IDX_test_2*']);
%% NEW head second part
Tail=[]; Head=[];
for i=1:size(files2,1)
   load([myKsDir,'\',files2(i).name]) 
   Head=[Head;HEAD]; Tail=[ Tail; TAIL];  
end

Tail(Head(:,1)>400,:)=[]; AUX=Head;
AUX(Head(:,1)>400,:)=[]; Head=AUX;


figure;
for i=1:size(Tail,1)
    line([Head(i,1);Tail(i,1)],[Head(i,2);Tail(i,2)],'LineWidth',1.75,'Color',[0 0.4470 0.7410 0.1] )
    hold on;
    plot(Head(i,1),Head(i,2),'.k','MarkerSize',10,'Color',[0 0 0 0.1]); axis equal
end

%%
Data=Head;
X = 0:1:752;
Y = 0:1:832;

D = length(Data(1,:));
Mu = mean(Data);
Sigma = cov(Data);
P_Gaussian = zeros(length(X),length(Y)); AUX=zeros(length(X),length(Y));
for i=1:length(X)
   for j1=1:length(Y)
       x = [X(i),Y(j1)];
       P_Gaussian(i,j1) = 1/((2*pi)^(D/2)*sqrt(det(Sigma)))...
                    *exp(-1/2*(x-Mu)*Sigma^-1*(x-Mu)');
     AUX(i,j1)=sum(Data(:,1)>=X(i) & Data(:,1)<X(i)+1 & Data(:,2)>=Y(j1) & Data(:,2)<Y(j1)+1);
   end
end
figure; contourb(P_Gaussian,100);
xlim([0 832]); ylim([0 758]); axis equal; colormap(brewermap(100,'GnBu'))

figure; contourb(AUX); 
xlim([0 832]); ylim([0 758]); axis equal; colormap(brewermap(100,'GnBu'))



%%
figure;
for i=1:size(Tail,1)
    line([Head(i,1);Tail(i,1)],[Head(i,2);Tail(i,2)],'LineWidth',1.75 )
    hold on;
    plot(Head(i,1),Head(i,2),'.k','MarkerSize',30)
end
axis equal;


ID=[]; SPE=[];
for i=1:size(files3,1)
   load([myKsDir,'\',files3(i).name]) 
   SPE=[SPE;FISHSP]; ID=[ ID; TAIL(:,3)];  
end

ID(SPE>300)=[]; SPE(SPE>300)=[]; 
figure; boxplot(SPE,ID,'PlotStyle','compact')


figure; hist3([mean([Head(:,1) Tail(:,1)],2) mean([Head(:,2) Tail(:,2)],2)],'Ctrs',{0:5:400 0:5:400},'CDataMode','auto','FaceColor','interp')
axis equal

figure; hist3(Head,[50 50],'CDataMode','auto','FaceColor','interp')
axis equal

