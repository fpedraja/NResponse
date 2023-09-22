%% head and tail run NEW IMPORTANT

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\8Fish_new_exp_2\');
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
myKsDir = uigetdir('Z:\locker\Fede\7Fish_new_exp\');
files2=dir([myKsDir, '\*_HEAD_IDX_control*']);
files3=dir([myKsDir, '\*_HEAD_IDX_test_1*']);
files4=dir([myKsDir, '\*_HEAD_IDX_test_2*']);
%% NEW head second part
Tail=[]; Head=[];
for i=1:size(files3,1)
   load([myKsDir,'\',files3(i).name]) 
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