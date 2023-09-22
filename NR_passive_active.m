clearvars; %close all; %addpath(genpath('C:\Users\fedu1\Downloads\sigTOOL'))
load(['Z:\locker\2FISH\20210415\20210415_007.mat']);
EODtime=V20210415_007_Ch2.values;  events=V20210415_007_Ch4.values;
%%
clearvars; %close all; %addpath(genpath('C:\Users\fedu1\Downloads\sigTOOL'))
load(['Z:\locker\2FISH\20210415\20210415_008.mat']);
EODtime=V20210415_008_Ch2.values;  events=V20210415_008_Ch4.values;
%%
clearvars; %close all; %addpath(genpath('C:\Users\fedu1\Downloads\sigTOOL'))
load(['Z:\locker\2FISH\20210415\20210415_006.mat']);
EODtime=V20210415_006_Ch2.values;  events=V20210415_006_Ch4.values;
%%
clearvars; %close all; %addpath(genpath('C:\Users\fedu1\Downloads\sigTOOL'))
load(['Z:\locker\2FISH\20210415\20210415_009.mat']);
EODtime=V20210415_009_Ch2.values;  events=V20210415_009_Ch4.values;
%%
clearvars; %close all; %addpath(genpath('C:\Users\fedu1\Downloads\sigTOOL'))
load(['Z:\locker\2FISH\20210415\20210415_010.mat']);
EODtime=V20210415_010_Ch2.values;  events=V20210415_010_Ch4.values;
%%
samplerate1=20000;
eventsamplerate=1000;

[value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.003,'MINPEAKDISTANCE',90);
[value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.04,'MINPEAKDISTANCE',200);
figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');

for k=1:size(sample4,1)
    sample3(sample3(:)==sample4(k))=[];
end


[value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',0.3,'MINPEAKDISTANCE',1900);
subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
sample1_N=sample1*(samplerate1/eventsamplerate);
EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1;
EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2;

% for i=1:size(sample1_N)
%     [idx2,idx]=min(abs(sample3-sample1_N(i)));
%     Freq1(i,:)=EODrate1(idx-10:idx+10);
% end
%  [xs, ys]=FitVal_EI([-10:10], mean(Freq1), [-10 10],0.9999);
%  [xs, ystd]=FitVal_EI([-10:10], mad(Freq1,1), [-10 10],0.9999);
%  figure;  addpath('D:\KIT3')
%  [hl, hp]=boundedline(xs', ys,ystd,'-k');

figure;a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
for i=1:size(sample1_N)
    [idx2,idx]=min(abs(sample3-sample1_N(i)));
    try
        Freq1post(a,:)=(sample3(idx-30:idx+30)-sample1_N(i))./samplerate1;
        Freq1(a,:)=EODr1(idx-30:idx+30);
        [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
        plot(xs,ys(a,:),'-k'); hold on;
        a=a+1;
    end
end
addpath('D:\KIT3'); addpath('C:\Users\fedu1\Google Drive\training')
% figure; [hl, hp]=boundedline(xs(:,1:end-1)', mean(diff(ys,1,2))',std(diff(ys,1,2))','-k');
figure; subplot(1,4,2);boxplot(EODr1,'PlotStyle','compact');
subplot(1,4,1); [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');

a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
for i=1:size(sample1_N)
    [idx2,idx]=min(abs(sample4-sample1_N(i)));
    try
        Freq1post(a,:)=(sample4(idx-90:idx+90)-sample1_N(i))./samplerate1;
        Freq1(a,:)=EODr2(idx-90:idx+90);
        [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
        a=a+1;
    end
end

% figure; [hl, hp]=boundedline(xs(:,1:end-1)', mean(diff(ys,1,2))',std(diff(ys,1,2))','-b');
subplot(1,4,4);boxplot(EODr2,'PlotStyle','compact');
subplot(1,4,3); [hl, hp]=boundedline(xs', mean(ys),std(ys),'-b');


%  [xs, ys]=FitVal_EI([-10:10], mean(Freq1), [-10 10],0.9999);
%  [xs, ystd]=FitVal_EI([-10:10], mad(Freq1,1), [-10 10],0.9999);
%  figure;  addpath('D:\KIT3')
%  [hl, hp]=boundedline(xs', ys,ystd,'-k');

%% control
clearvars; %addpath(genpath('C:\Users\fedu1\Downloads\sigTOOL'))
load(['Z:\locker\2FISH\20210415\20210415_011.mat']);
EODtime=V20210415_011_Ch2.values;  events=V20210415_011_Ch4.values; MINVAL1=[]; MINVAL2=[];
%%
clearvars; %addpath(genpath('C:\Users\fedu1\Downloads\sigTOOL'))
load(['Z:\locker\2FISH\20210415\20210415_004.mat']);
EODtime=V20210415_004_Ch2.values;  events=V20210415_004_Ch4.values; MINVAL1=[]; MINVAL2=[];
%%
clearvars; %addpath(genpath('C:\Users\fedu1\Downloads\sigTOOL'))
load(['Z:\locker\2FISH\20210415\20210415_012.mat']);
EODtime=V20210415_012_Ch2.values;  events=V20210415_012_Ch4.values; MINVAL1=[]; MINVAL2=[];
%%
samplerate1=20000;
eventsamplerate=1000;

[value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.003,'MINPEAKDISTANCE',90);
[value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',0.3,'MINPEAKDISTANCE',1900);
sample1_N=sample1*(samplerate1/eventsamplerate);
EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1;
figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');


a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
for i=1:size(sample1_N)
    [idx2,idx]=min(abs(sample3-sample1_N(i)));
    try
        Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
        Freq1(a,:)=EODr1(idx-40:idx+40);
        [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
        a=a+1;
    end
    %plot(Freq1(i,:),i,'.k'); hold on;
    
end
figure; subplot(1,4,2);boxplot(EODr1,'PlotStyle','compact');
subplot(1,4,1); [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');

a=1;
for g=180:10:340
    AUX=[]; AUX1=[];
    AUX=reshape(ys2(:,g:g+9),[],1);
    AUX1=reshape(ys(:,g:g+9),[],1);
    GR=[zeros(size(AUX,1),1); ones(size(AUX1,1),1)];
    panova(a) = anova1([AUX;AUX1],GR,'off');  tim(a)=xs(1,g);
    a=a+1;
end
figure; plot(tim,panova)

%% NR after second fish EOD after OBJ on

MINVAL1=[]; MINVAL2=[];
samplerate1=20000;
eventsamplerate=1000;

[value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.003,'MINPEAKDISTANCE',90);
[value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.04,'MINPEAKDISTANCE',200);
figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');

for k=1:size(sample4,1)
    sample3(sample3(:)==sample4(k))=[];
end


[value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',0.3,'MINPEAKDISTANCE',1900);
subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
sample1_N=sample1*(samplerate1/eventsamplerate);
EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1;
EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2;

% for i=1:size(sample1_N)
%     [idx2,idx]=min(abs(sample3-sample1_N(i)));
%     Freq1(i,:)=EODrate1(idx-10:idx+10);
% end
%  [xs, ys]=FitVal_EI([-10:10], mean(Freq1), [-10 10],0.9999);
%  [xs, ystd]=FitVal_EI([-10:10], mad(Freq1,1), [-10 10],0.9999);
%  figure;  addpath('D:\KIT3')
%  [hl, hp]=boundedline(xs', ys,ystd,'-k');

a=1; xs=[]; ys=[]; xt=[]; yt=[];
for i=1:size(sample1_N)
    [idx2,idx]=min(abs(sample4-sample1_N(i)));
    if sample4(idx)>sample1_N(i)
        [idx3,idx4]=min(abs(sample3-sample4(idx)));
        try
            Freq1post(a,:)=(sample3(idx4-50:idx4+50)-sample4(idx))./samplerate1;
            Freq1(a,:)=EODr1(idx4-50:idx4+50);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            a=a+1;
        end
    else
        [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
        try
            Freq1post(a,:)=(sample3(idx4-50:idx4+50)-sample4(idx))./samplerate1;
            Freq1(a,:)=EODr1(idx4-50:idx4+50);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            a=a+1;
        end
    end
    %plot(Freq1(i,:),i,'.k'); hold on;
    
    
    
end
addpath('D:\KIT3'); addpath('C:\Users\fedu1\Google Drive\training')
% figure; [hl, hp]=boundedline(xs(:,1:end-1)', mean(diff(ys,1,2))',std(diff(ys,1,2))','-k');
figure; subplot(1,2,1); [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');
%% random bootsrap
samplerate1=20000;
eventsamplerate=1000;

[value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.003,'MINPEAKDISTANCE',90);
[value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.04,'MINPEAKDISTANCE',200);
figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');

for k=1:size(sample4,1)
    sample3(sample3(:)==sample4(k))=[];
end


[value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',0.3,'MINPEAKDISTANCE',1900);
subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
sample1_N=sample1*(samplerate1/eventsamplerate);
EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1;
EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2;

% for i=1:size(sample1_N)
%     [idx2,idx]=min(abs(sample3-sample1_N(i)));
%     Freq1(i,:)=EODrate1(idx-10:idx+10);
% end
%  [xs, ys]=FitVal_EI([-10:10], mean(Freq1), [-10 10],0.9999);
%  [xs, ystd]=FitVal_EI([-10:10], mad(Freq1,1), [-10 10],0.9999);
%  figure;  addpath('D:\KIT3')
%  [hl, hp]=boundedline(xs', ys,ystd,'-k');

b=1;
for j=1:100000
    xs=[]; ys=[]; Freq1post=[]; Freq1=[]; a=1;
    sample1_N=randi(max(sample3),75,1);
    for i=1:75
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            %plot(xs,ys(a,:),'-k'); hold on;
            a=a+1;
        end
    end
    try
        pkruskal(b) = kruskalwallis([AUX(1:73) ys(1:73,268)],[],'off');
        panova(b) = anova1([AUX(1:73) ys(1:73,268)],[],'off');
        b=b+1;
    end
end
addpath('D:\KIT3'); addpath('C:\Users\fedu1\Google Drive\training')
% figure; [hl, hp]=boundedline(xs(:,1:end-1)', mean(diff(ys,1,2))',std(diff(ys,1,2))','-k');
subplot(1,4,1); [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');
%% bootstrapping 2
samplerate1=20000;
eventsamplerate=1000;

[value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.003,'MINPEAKDISTANCE',90);
[value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.04,'MINPEAKDISTANCE',200);
figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');

for k=1:size(sample4,1)
    sample3(sample3(:)==sample4(k))=[];
end


[value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',0.3,'MINPEAKDISTANCE',1900);
subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
sample1_N=sample1*(samplerate1/eventsamplerate);
EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1;
EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2;



for j=1:10000
    xs=[]; ys=[]; Freq1post=[]; Freq1=[]; a=1;
    sample1_N=randi(max(sample3),75,1);
    for i=1:75
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-30:idx+30)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-30:idx+30);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            %plot(xs,ys(a,:),'-k'); hold on;
            b=1;
            for g=180:10:340
                AUX=[]; AUX1=[];
                AUX=reshape(ys2(:,g:g+9),[],1);
                AUX1=reshape(ys(:,g:g+9),[],1);
                GR=[zeros(size(AUX,1),1); ones(size(AUX1,1),1)];
                panova(a,b) = anova1([AUX;AUX1],GR,'off');  tim(b)=xs(1,g);
                b=b+1;
            end
            a=a+1;
        end
    end
    %     try
    %             pkruskal(b) = kruskalwallis([AUX(1:73) ys(1:73,268)],[],'off');
    %             panova(b) = anova1([AUX(1:73) ys(1:73,268)],[],'off');
    %             b=b+1;
    %     end
end
figure; plot(tim,mean(panova))
figure; plot(tim,median(panova))
%% NR calculation with z-score

samplerate1=20000;
eventsamplerate=1000;

[value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.003,'MINPEAKDISTANCE',90);
[value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.04,'MINPEAKDISTANCE',200);
figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');

for k=1:size(sample4,1)
    sample3(sample3(:)==sample4(k))=[];
end


[value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',0.3,'MINPEAKDISTANCE',1900);
subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
sample1_N=sample1*(samplerate1/eventsamplerate);
EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=diff(EODr1); EODr1(end+1)=EODr1(end); EODr1=zscore(EODr1);
EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2=diff(EODr2); EODr2(end+1)=EODr2(end); EODr2=zscore(EODr2);

% for i=1:size(sample1_N)
%     [idx2,idx]=min(abs(sample3-sample1_N(i)));
%     Freq1(i,:)=EODrate1(idx-10:idx+10);
% end
%  [xs, ys]=FitVal_EI([-10:10], mean(Freq1), [-10 10],0.9999);
%  [xs, ystd]=FitVal_EI([-10:10], mad(Freq1,1), [-10 10],0.9999);
%  figure;  addpath('D:\KIT3')
%  [hl, hp]=boundedline(xs', ys,ystd,'-k');

figure;a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
for i=1:size(sample1_N)
    [idx2,idx]=min(abs(sample3-sample1_N(i)));
    try
        Freq1post(a,:)=(sample3(idx-30:idx+30)-sample1_N(i))./samplerate1;
        Freq1(a,:)=EODr1(idx-30:idx+30);
        [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
        plot(xs,ys(a,:),'-k'); hold on;
        a=a+1;
    end
end
addpath('D:\KIT3'); addpath('C:\Users\fedu1\Google Drive\training')
% figure; [hl, hp]=boundedline(xs(:,1:end-1)', mean(diff(ys,1,2))',std(diff(ys,1,2))','-k');
figure; subplot(1,4,2);boxplot(EODr1,'PlotStyle','compact');
subplot(1,4,1); [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');

a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
for i=1:size(sample1_N)
    [idx2,idx]=min(abs(sample4-sample1_N(i)));
    try
        Freq1post(a,:)=(sample4(idx-90:idx+90)-sample1_N(i))./samplerate1;
        Freq1(a,:)=EODr2(idx-90:idx+90);
        [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
        a=a+1;
    end
end

% figure; [hl, hp]=boundedline(xs(:,1:end-1)', mean(diff(ys,1,2))',std(diff(ys,1,2))','-b');
subplot(1,4,4);boxplot(EODr2,'PlotStyle','compact');
subplot(1,4,3); [hl, hp]=boundedline(xs', mean(ys),std(ys),'-b');
%% 23-04-21
%NEW DATA NR 2 FISH

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\Fede\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% first part, 1 fish

objdist=[10 0 5 60  25 31.6  25 31.6];
vidnum=[1 2 3 4 6 7 19 20];



samplerate1=20000;
eventsamplerate=1000;
%%
figure;
for j=1:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1;   %EODr1=zscore(EODr1); %EODr1(1:end-1)=EODr1(2:end);   EODr1(end+1)=EODr1(end);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    %figure; subplot 131
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
           % plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %subplot 133; boxplot(EODr1,'PlotStyle','compact');
    subplot (1,length(vidnum),j); [hl, hp]=boundedline(xs', mean(ys),std(ys)/sqrt(size(ys,1)),'-k'); xlim([-1 2.5]); ylim([0 20]);
    
    freqDIFF(j)=(mean(mean(ys(:,176:201))))-(mean(mean(ys(:,126:151))));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
end


% figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 27 -0.2]);
% figure; [param, stat]=sigm_fit(objdist,perf,[nan nan 27 -0.2]);
% R1=1-(sum((stat.ypred-perf').^2)/sum((mean(perf)-perf).^2));
%% %% NEW DATA NR 2 FISH second part

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\Fede\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

objdist=[31.6 60 40 0 25 10];
vidnum=[9 10 11 12 13 14];

samplerate1=20000;
eventsamplerate=1000;
figure;
for j=1:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    sample3=[]; value3=[]; sample4=[]; value4=[];
    [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    a=1; b=1;
    for k=1:length(value2)
        if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
            sample3(a)=sample2(k);
            value3(a)=value2(k);
            a=a+1;
        else
            sample4(b)=sample2(k);
            value4(b)=value2(k);
            b=b+1;
        end
    end
    
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1;  EODr1(2:end+1)=EODr1; %EODr1=zscore(EODr1);
    % figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    % subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[]; ys1=[];
    %figure; subplot 141
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            [xs(a,:), ys1(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
           % plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            if NRtrialsfish2(i)==0
                % plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            end
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    %subplot 142; 
     subplot (1,length(vidnum),j); [hl, hp]=boundedline(xs', mean(ys),std(ys)/sqrt(size(ys,1)),'-k'); xlim([-1 2.5]); ylim([0 20]);
    freqDIFF(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    NRtrialsfish1=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    
% second fish
    %     subplot 133; boxplot(EODr1,'PlotStyle','compact');
    %      subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2(2:end+1)=EODr2; %EODr2=zscore(EODr2);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[]; ys2=[];
    %subplot 143
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        try
            Freq1post(a,:)=(sample4(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr2(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            [xs(a,:), ys2(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            %plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    freqDIFFfish2(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perffish2(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
   % subplot 144; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    NRtrialsfish2=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    
    a=1; Interv=[];
    for i=1:size(sample3,2)
        [idx2,idx]=min(abs(sample4-sample3(i)));
        
        if sample4(idx)-sample3(i)<0
            idx=idx+1;
        end
        
        try
            Interv(a)=(sample4(idx)-sample3(i))/samplerate1;
        catch
            Interv(a)=0;
            %plot(Freq1(a-1,:),i,'.k'); hold on;
        end
        if Interv(a)<=0
            Interv(a)=0;
        end
        a=a+1;
    end
    
    %
    %figure; violinplot(Interv)
    %figure; boxplot(Interv,'PlotStyle','compact')
    %     figure; hist(Interv,500)
    %
    %     figure; plot(sample3/samplerate1,Interv,'.k')
   % hold on;
    %plot((1:1:size(events,1))/eventsamplerate,events);
    
    AUX=[]; AUX2=sample3(Interv>=0.10 & Interv<=0.20)/samplerate1; AUX3=[];
    for l=1:size(sample1_N,1)
        AUX3(l)=sum(AUX2>=sample1_N(l)/samplerate1 & AUX2<=(sample1_N(l)/samplerate1)+1);
    end
    AUX3(AUX3==0)=[];
    perfInterv(j)=size(AUX3,2)/size(sample1_N,1);
end


% figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 31 -0.2]);
% figure; [param,stat1]=sigm_fit(objdist,perf,[nan nan 33 -0.2]);
% R2_1=1-(sum((stat1.ypred-perf').^2)/sum((mean(perf)-perf).^2));
%
% figure; [param]=sigm_fit(objdist,freqDIFFfish2,[NaN NaN 31 -0.2]);
% figure; [param,stat2]=sigm_fit(abs(objdist),perffish2,[nan nan NaN -0.2]);
% R2_2=1-(sum((stat2.ypred-perffish2').^2)/sum((mean(perffish2)-perffish2).^2));
% % figure; plot(objdist,freqDIFFfish2,'.k');
% [param]=sigm_fit(objdist,freqDIFFfish2,[NaN NaN NaN NaN]);


%figure; [param,stat1]=sigm_fit(objdist,perfInterv,[nan nan 22 -0.2]);
%R2_Interv=1-(sum((stat1.ypred-perfInterv').^2)/sum((mean(perfInterv)-perfInterv).^2));
%% %%%%%%%% NR after EODs of second fish (after obj on)

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% % first part, 1 fish
% objdist=[0 20 22.36 60.8 51 41.23 31.6 10 26.9 36.4 0 23.25];
% vidnum=[1 2 3 4 5 6 7 8 9 10 11 12];
objdist=[31.6 60 40 0 25 10];
vidnum=[9 10 11 12 13 14];

samplerate1=20000;
eventsamplerate=1000;
%%
for j=1%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    sample3=[]; value3=[]; sample4=[]; value4=[];
    [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.003,'MINPEAKDISTANCE',90);
    a=1; b=1;
    for k=1:length(value2)
        if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
            sample3(a)=sample2(k);
            value3(a)=value2(k);
            a=a+1;
        else
            sample4(b)=sample2(k);
            value4(b)=value2(k);
            b=b+1;
        end
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2=zscore(EODr2);
    
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+3)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+3))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
end


%% 30-04-21
%%%% NEW DATA NR 2 FISH

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% first part, 1 fish
objdist=[20 22.36 60.8 51 41.23 31.6 26.9 36.4 0 23.25];
vidnum=[2 3 4 5 6 7 9 10 11 12];



samplerate1=20000;
eventsamplerate=1000;
%%
for j=1%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0002,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1); %EODr1(1:end-1)=EODr1(2:end);   EODr1(end+1)=EODr1(end);
    %     figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    %     subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    %
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 131
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999999999);
            plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    subplot 133; boxplot(EODr1,'PlotStyle','compact');
    subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    
    freqDIFF(j)=(mean(mean(ys(:,176:201))))-(mean(mean(ys(:,126:151))));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
end


figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 27 -0.2]);
figure; [param, stat]=sigm_fit(objdist,perf,[0.55 nan nan -0.2]);
R1=1-(sum((stat.ypred-perf').^2)/sum((mean(perf)-perf).^2));

%% second part


clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23  30];
vidnum=[13 14 15 16 17 18 19 20 21 22  25];

samplerate1=20000;
eventsamplerate=1000;
perfInterv=[];
for j=4%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0006,'MINPEAKDISTANCE',200);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 141
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999999);
            plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
            plot(Freq1(a-1,:),i,'.k'); hold on;
        end
        
        
    end
    freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    subplot 142; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    freqDIFF(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    NRtrialsfish2=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    
    % second fish
    %     subplot 133; boxplot(EODr1,'PlotStyle','compact');
    %      subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2=zscore(EODr2);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    % subplot 143
    %figure;
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        try
            Freq1post(a,:)=(sample4(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr2(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            %             if NRtrialsfish2(i)==0
            %             plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            %             end
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    freqDIFFfish2(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perffish2(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    %     subplot 144; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    NRtrialsfish1=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    
    a=1; Interv=[];
    for i=1:size(sample4,1)
        [idx2,idx]=min(abs(sample3-sample4(i)));
        
        if sample3(idx)-sample4(i)<0
            idx=idx+1;
        end
        
        try
            Interv(a)=(sample3(idx)-sample4(i))/samplerate1;
        catch
            Interv(a)=0;
            %plot(Freq1(a-1,:),i,'.k'); hold on;
        end
        if Interv(a)<=0
            Interv(a)=0;
        end
        a=a+1;
    end
    %figure; hist(Interv,500)
    AUX=[]; AUX2=[]; AUX2=sample4(Interv>=0.013 & Interv<=0.020)/samplerate1; AUX3=[];
    for l=1:size(sample1_N,1)
        AUX3(l)=sum(AUX2>=sample1_N(l)/samplerate1 & AUX2<=(sample1_N(l)/samplerate1)+1);
    end
    AUX3(AUX3==0)=[];
    perfInterv(j)=size(AUX3,2)/size(sample1_N,1);
    
    
end


% figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 31 -0.2]);
% figure; [param,stat1]=sigm_fit(objdist,perf,[nan nan nan -0.2]);
% % R2_1=1-(sum((stat1.ypred-perf').^2)/sum((mean(perf)-perf).^2));
% perffish2=perffish2/0.9;
% figure; [param]=sigm_fit(objdist,freqDIFFfish2,[NaN NaN 31 -0.2]);
% figure; [param,stat2]=sigm_fit(abs(objdist),perffish2,[nan nan 28 -0.1]);
% R2_2=1-(sum((stat2.ypred-perffish2').^2)/sum((mean(perffish2)-perffish2).^2));

figure; [param,stat1]=sigm_fit(objdist,perfInterv,[nan nan nan -0.5]);
R2_Interv=1-(sum((stat1.ypred-perfInterv').^2)/sum((mean(perfInterv)-perfInterv).^2));
%% test
j=4;

load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0006,'MINPEAKDISTANCE',200);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2(2:end+1)=EODr2; EODr2=zscore(EODr2);
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+3)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+3))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])


%% 07-05-2021
%%%% NEW DATA NR 2 FISH

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% first part, 1 fish
%objdist=[20 22.36 60.8 51 41.23 31.6 26.9 36.4 0 23.25];
objdist=[0 5 51 31.6 22.36 10 15 41.23 60.8 20 26.9 23.25 0];
vidnum=[1 2 3 4 5 6 7 8 9 10 11 12 13];



samplerate1=20000;
eventsamplerate=1000;

for j=2%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0002,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1); %EODr1(1:end-1)=EODr1(2:end);   EODr1(end+1)=EODr1(end);
    %  figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    % subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    %
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 131
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999999999999);
            plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %9         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    subplot 133; boxplot(EODr1,'PlotStyle','compact');
    subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    
    freqDIFF(j)=(mean(mean(ys(:,176:201))))-(mean(mean(ys(:,126:151))));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
end


figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 27 -0.2]);
figure; [param, stat]=sigm_fit(objdist,perf,[nan nan nan -0.2]);
R1=1-(sum((stat.ypred-perf').^2)/sum((mean(perf)-perf).^2));

%% second part


clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

%objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23 10 30];
objdist=[0 51 22.36 15 0 20 26.9 60.8 41.23  22.36 23.25 31.6 5];
vidnum=[14 15 16 17 18 19 20 21 22 23 24 25 26];

samplerate1=20000;
eventsamplerate=1000;
perfInterv=[];
for j=6%1:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.002,'MINPEAKDISTANCE',200);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
    figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 141
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999999999);
            plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    subplot 142; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    freqDIFF(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    NRtrialsfish2=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    
    % second fish
    %     subplot 133; boxplot(EODr1,'PlotStyle','compact');
    %      subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2=zscore(EODr2);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    % subplot 143
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        try
            Freq1post(a,:)=(sample4(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr2(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            % plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    freqDIFFfish2(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perffish2(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    %     subplot 144; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    NRtrialsfish1=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    a=1; Interv=[];
    for i=1:size(sample4,1)
        [idx2,idx]=min(abs(sample3-sample4(i)));
        
        if sample3(idx)-sample4(i)<0
            idx=idx+1;
        end
        
        try
            Interv(a)=(sample3(idx)-sample4(i))/samplerate1;
        catch
            Interv(a)=0;
            %plot(Freq1(a-1,:),i,'.k'); hold on;
        end
        if Interv(a)<=0
            Interv(a)=0;
        end
        a=a+1;
    end
    %figure; hist(Interv,500)
    AUX=[]; AUX2=[]; AUX2=sample4(Interv>=0.013 & Interv<=0.020)/samplerate1; AUX3=[];
    for l=1:size(sample1_N,1)
        AUX3(l)=sum(AUX2>=sample1_N(l)/samplerate1 & AUX2<=(sample1_N(l)/samplerate1)+1);
    end
    AUX3(AUX3==0)=[];
    perfInterv(j)=size(AUX3,2)/size(sample1_N,1);
    
end

%  figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 31 -0.2]);
%  figure; [param,stat1]=sigm_fit(objdist,perf,[nan nan nan -0.2]);
% % R2_1=1-(sum((stat1.ypred-perf').^2)/sum((mean(perf)-perf).^2));
% perffish2=perffish2*0.9;
% figure; [param]=sigm_fit(objdist,freqDIFFfish2,[NaN NaN 31 -0.2]);
% figure; [param,stat2]=sigm_fit(abs(objdist),perffish2,[nan nan nan nan]);
% R2_2=1-(sum((stat2.ypred-perffish2').^2)/sum((mean(perffish2)-perffish2).^2));
figure; [param,stat1]=sigm_fit(objdist,perfInterv,[nan nan nan -0.6]);
R2_Interv=1-(sum((stat1.ypred-perfInterv').^2)/sum((mean(perfInterv)-perfInterv).^2));
%% test 07-05-2021
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

%objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23 10 30];
objdist=[0 51 22.36 15 0 20 26.9 60.8 41.23  22.36 23.25 31.6 5];
vidnum=[14 15 16 17 18 19 20 21 22 23 24 25 26];

samplerate1=20000;
eventsamplerate=1000;
j=6;%1:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.002,'MINPEAKDISTANCE',200);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
     EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2;  EODr2(2:end+1)=EODr2; EODr2=zscore(EODr2);
    
    

    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+3)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+3))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])

%% 13-05-2021
%%%% NEW DATA NR 2 FISH

% first part, 1 fish
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);


%objdist=[20 22.36 60.8 51 41.23 31.6 26.9 36.4 0 23.25];
objdist=[0 60.8 41.23 51 31.6 5 15 20 10 26.9 0 23.25 10 12 7];
vidnum=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];



samplerate1=20000;
eventsamplerate=1000;
%%
for j=14%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0002,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1); %EODr1(1:end-1)=EODr1(2:end);   EODr1(end+1)=EODr1(end);
    figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    % subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    %
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 131
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.99999999999);
            plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    subplot 133; boxplot(EODr1,'PlotStyle','compact');
    subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    
    freqDIFF(j)=(mean(mean(ys(:,176:201))))-(mean(mean(ys(:,126:151))));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
end


figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 27 -0.2]);
figure; [param, stat]=sigm_fit(objdist,perf,[nan nan nan -0.3]);
R1=1-(sum((stat.ypred-perf').^2)/sum((mean(perf)-perf).^2));

%% second part

close all
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

%objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23 10 30];
objdist=[7 0 20 31.6 41.23 7 10 5 23.25  12 51 60.8];
vidnum=[16 17 19 20 21 22 23 24 25 27 28 29];

samplerate1=20000;
eventsamplerate=1000;
%%
for j=3%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    if j>=3 && j<=4
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.01,'MINPEAKDISTANCE',200);
    elseif j==7
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.009,'MINPEAKDISTANCE',200);
    elseif j==12
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.013,'MINPEAKDISTANCE',200);
    else
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.008,'MINPEAKDISTANCE',200);
    end
    figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 141
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.99999999999999);
            plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    subplot 142; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    freqDIFF(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    
    
    % second fish
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     subplot 133; boxplot(EODr1,'PlotStyle','compact');
    %      subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2=zscore(EODr2);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    % subplot 143
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        try
            Freq1post(a,:)=(sample4(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr2(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            % plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    freqDIFFfish2(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perffish2(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    %     subplot 144; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    
    
end


figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 31 -0.2]);
figure; [param,stat1]=sigm_fit(objdist,perf,[nan nan nan -0.2]);
% R2_1=1-(sum((stat1.ypred-perf').^2)/sum((mean(perf)-perf).^2));
figure; [param]=sigm_fit(objdist,freqDIFFfish2,[NaN NaN 31 -0.2]);

% perffish2(perffish2>=0.35)=perffish2(perffish2>=0.35)/0.7;
figure; [param,stat2]=sigm_fit(objdist,perffish2,[nan nan nan -0.2]);
R2_2=1-(sum((stat2.ypred-perffish2').^2)/sum((mean(perffish2)-perffish2).^2));


%% test 13-05-2021
close all
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

%objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23 10 30];
objdist=[7 0 20 31.6 41.23 7 10 5 23.25  12 51 60.8];
vidnum=[16 17 19 20 21 22 23 24 25 27 28 29];

samplerate1=20000;
eventsamplerate=1000;
j=3;%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    if j>=3 && j<=4
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.01,'MINPEAKDISTANCE',200);
    elseif j==7
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.009,'MINPEAKDISTANCE',200);
    elseif j==12
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.013,'MINPEAKDISTANCE',200);
    else
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.008,'MINPEAKDISTANCE',200);
    end
   % figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2;  EODr2(2:end+1)=EODr2; EODr2=zscore(EODr2);
    
    

    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+3)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+3))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
%% 12-05-2021
%%%% NEW DATA NR 2 FISH

% first part, 1 fish
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);


%objdist=[20 22.36 60.8 51 41.23 31.6 26.9 36.4 0 23.25];
objdist=[0 5 51 60.8 31.6 10 22.36 31.6 41.23 15 26.9 0];
vidnum=[1 2 3 4 5 6 7 8 9 10 11 12];



samplerate1=20000;
eventsamplerate=1000;

for j=1:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0002,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1); %EODr1(1:end-1)=EODr1(2:end);   EODr1(end+1)=EODr1(end);
    figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    % subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    %
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    %     figure; subplot 131
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            % plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %subplot 133; boxplot(EODr1,'PlotStyle','compact');
    %     subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    
    freqDIFF(j)=(mean(mean(ys(:,176:201))))-(mean(mean(ys(:,126:151))));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
end


figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 27 -0.2]);
figure; [param, stat]=sigm_fit(objdist,perf,[nan nan 5 -0.2]);
R1=1-(sum((stat.ypred-perf').^2)/sum((mean(perf)-perf).^2));

%% second part


clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

%objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23 10 30];
objdist=[0 31.6 41.23 51 60.8 20 26.9 15 5 15 22.36];
vidnum=[13 14 15 16 17 18 19 20 21 25 26];

samplerate1=20000;
eventsamplerate=1000;
perfInterv=[];
for j=1%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    if j>=3 && j<=4
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.01,'MINPEAKDISTANCE',200);
    else
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.008,'MINPEAKDISTANCE',200);
    end
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    %     figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %     subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    %
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    % figure; subplot 141
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            %plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    % subplot 142; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    freqDIFF(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    NRtrialsfish2=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    
    % second fish
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     subplot 133; boxplot(EODr1,'PlotStyle','compact');
    %      subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2=zscore(EODr2);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    % subplot 143
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        try
            Freq1post(a,:)=(sample4(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr2(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            % plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    freqDIFFfish2(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perffish2(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    %     subplot 144; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    NRtrialsfish1=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    
    a=1; Interv=[];
    for i=1:size(sample4,1)
        [idx2,idx]=min(abs(sample3-sample4(i)));
        
        if sample3(idx)-sample4(i)<0
            idx=idx+1;
        end
        
        try
            Interv(a)=(sample3(idx)-sample4(i))/samplerate1;
        catch
            Interv(a)=0;
            %plot(Freq1(a-1,:),i,'.k'); hold on;
        end
        if Interv(a)<=0
            Interv(a)=0;
        end
        a=a+1;
    end
    %figure; hist(Interv,500)
    AUX=[]; AUX2=[]; AUX2=sample4(Interv>=0.010 & Interv<=0.020)/samplerate1; AUX3=[];
    for l=1:size(sample1_N,1)
        AUX3(l)=sum(AUX2>=sample1_N(l)/samplerate1 & AUX2<=(sample1_N(l)/samplerate1)+0.5);
    end
    AUX3(AUX3==0)=[];
    perfInterv(j)=size(AUX3,2)/size(sample1_N,1);
    
end


%  figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 31 -0.2]);
%  figure; [param,stat1]=sigm_fit(objdist,perf,[nan nan nan -0.2]);
% % R2_1=1-(sum((stat1.ypred-perf').^2)/sum((mean(perf)-perf).^2));
% figure; [param]=sigm_fit(objdist,freqDIFFfish2,[NaN NaN 31 -0.2]);
%
% figure; [param,stat2]=sigm_fit(objdist,perffish2,[nan nan nan -0.3]);
% R2_2=1-(sum((stat2.ypred-perffish2').^2)/sum((mean(perffish2)-perffish2).^2));

figure; [param,stat1]=sigm_fit(objdist,perfInterv,[nan nan nan -0.3]);
R2_Interv=1-(sum((stat1.ypred-perfInterv').^2)/sum((mean(perfInterv)-perfInterv).^2));

%% test 




%% 21-05-2021
%%%% NEW DATA NR 2 FISH

% first part, 1 fish
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);


%objdist=[20 22.36 60.8 51 41.23 31.6 26.9 36.4 0 23.25];
objdist=[0 60.8 51 31.6 20 26.9 10 41.23 15 5 23.25 0];
vidnum=[1 2 3 4 5 6 7 8 9 10 11 12];



samplerate1=20000;
eventsamplerate=1000;
%%
for j=12%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0002,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1); %EODr1(1:end-1)=EODr1(2:end);   EODr1(end+1)=EODr1(end);
    figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    % subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    %
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 131
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999999999);
            plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    subplot 133; boxplot(EODr1,'PlotStyle','compact');
    subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    
    freqDIFF(j)=(mean(mean(ys(:,176:201))))-(mean(mean(ys(:,126:151))));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
end


figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN nan -0.2]);
figure; [param, stat]=sigm_fit(objdist,perf,[nan nan nan -0.15]);
R1=1-(sum((stat.ypred-perf').^2)/sum((mean(perf)-perf).^2));

%% second part


clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

%objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23 10 30];
objdist=[0 31.6  51 41.23 15 10 20 23.25 5 26.9];
vidnum=[13 15  18 19 20 21 22 23 25 26];

samplerate1=20000;
eventsamplerate=1000;
perfInterv=[];
%%
for j=7%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    if j>=3 && j<=4
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.01,'MINPEAKDISTANCE',200);
    else
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.008,'MINPEAKDISTANCE',200);
    end
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 121
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999999999);
            %plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    subplot 142; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    freqDIFF(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    NRtrialsfish2=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    
    % second fish
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     subplot 133; boxplot(EODr1,'PlotStyle','compact');
    %      subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2=zscore(EODr2);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    % subplot 122
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        try
            Freq1post(a,:)=(sample4(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr2(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            %plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    freqDIFFfish2(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perffish2(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
         subplot 144; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    NRtrialsfish1=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    a=1; Interv=[];
    for i=1:size(sample4,1)
        [idx2,idx]=min(abs(sample3-sample4(i)));
        
        if sample3(idx)-sample4(i)<0
            idx=idx+1;
        end
        
        try
            Interv(a)=(sample3(idx)-sample4(i))/samplerate1;
        catch
            Interv(a)=0;
            %plot(Freq1(a-1,:),i,'.k'); hold on;
        end
        if Interv(a)<=0
            Interv(a)=0;
        end
        a=a+1;
    end
    %figure; hist(Interv,500)
    AUX=[]; AUX2=[]; AUX2=sample4(Interv>=0.012 & Interv<=0.018)/samplerate1; AUX3=[];
    for l=1:size(sample1_N,1)
        AUX3(l)=sum(AUX2>=sample1_N(l)/samplerate1 & AUX2<=(sample1_N(l)/samplerate1)+0.5);
    end
    AUX3(AUX3==0)=[];
    perfInterv(j)=size(AUX3,2)/size(sample1_N,1);
    
end


%  figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 31 -0.2]);
%  figure; [param,stat1]=sigm_fit(objdist,perf,[nan nan nan -0.2]);
% % R2_1=1-(sum((stat1.ypred-perf').^2)/sum((mean(perf)-perf).^2));
% figure; [param]=sigm_fit(objdist,freqDIFFfish2,[NaN NaN 31 -0.2]);
%
%
% figure; [param,stat2]=sigm_fit(objdist,perffish2,[nan nan nan -0.2]);
% R2_2=1-(sum((stat2.ypred-perffish2').^2)/sum((mean(perffish2)-perffish2).^2));

figure; [param,stat1]=sigm_fit(objdist,perfInterv,[nan nan 28 -0.1]);
R2_Interv=1-(sum((stat1.ypred-perfInterv').^2)/sum((mean(perfInterv)-perfInterv).^2));

%% test 21-05-2021

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\Fede\2FISH\20210521\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

%objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23 10 30];
objdist=[0 31.6  51 41.23 15 10 20 23.25 5 26.9];
vidnum=[13 15  18 19 20 21 22 23 25 26];

samplerate1=20000;
eventsamplerate=1000;
perfInterv=[];

j=7;
load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    if j>=3 && j<=4
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.01,'MINPEAKDISTANCE',200);
    else
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.008,'MINPEAKDISTANCE',200);
    end
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
        EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2;  EODr2(2:end+1)=EODr2; EODr2=zscore(EODr2);
    
    

    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+3)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+3))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
%% 20-05-2021
%%%% NEW DATA NR 2 FISH

% first part, 1 fish
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);


%objdist=[20 22.36 60.8 51 41.23 31.6 26.9 36.4 0 23.25];
objdist=[0 60.8 51 41.23 26.9 15 10 31.6  20];
vidnum=[1 2 3 4 5 6 7 8 11 ];



samplerate1=20000;
eventsamplerate=1000;

for j=8%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0002,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1); %EODr1(1:end-1)=EODr1(2:end);   EODr1(end+1)=EODr1(end);
    % figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    % subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    %
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 131
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999999999);
            plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    subplot 133; boxplot(EODr1,'PlotStyle','compact');
    subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    
    freqDIFF(j)=(mean(mean(ys(:,176:201))))-(mean(mean(ys(:,126:151))));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.9)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
end

figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN nan -0.2]);
figure; [param, stat]=sigm_fit(objdist,perf,[nan nan 15 -0.1]);
R1=1-(sum((stat.ypred-perf').^2)/sum((mean(perf)-perf).^2));

%% second part


clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

%objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23 10 30];
objdist=[60.8 20 31.6 41.23 15 10 51 26.9];
vidnum=[15 17 18 19 20 21 22 25];

samplerate1=20000;
eventsamplerate=1000;
perfInterv=[];
for j=3%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    %     sample3=[]; value3=[]; sample4=[]; value4=[];
    %     [value2,sample2]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.03,'MINPEAKDISTANCE',90);
    %     a=1; b=1;
    %     for k=1:length(value2)
    %         if EODtime(sample2(k)+6)-EODtime(sample2(k)+19)<=-0.0015
    %             sample3(a)=sample2(k);
    %             value3(a)=value2(k);
    %             a=a+1;
    %         else
    %             sample4(b)=sample2(k);
    %             value4(b)=value2(k);
    %             b=b+1;
    %         end
    %     end
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',90);
    if j>=3 && j<=4
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.005,'MINPEAKDISTANCE',200);
    else
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.005,'MINPEAKDISTANCE',200);
    end
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
    % figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    % subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; subplot 141
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            %plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    subplot 142; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    freqDIFF(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perf(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    NRtrialsfish2=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    
    
    % second fish
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     subplot 133; boxplot(EODr1,'PlotStyle','compact');
    %      subplot 132; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2; EODr2=zscore(EODr2);
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    %subplot(2,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    
    
    a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    % subplot 143
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        try
            Freq1post(a,:)=(sample4(idx-40:idx+40)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr2(idx-40:idx+40);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
            % plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
            a=a+1;
        end
        %         plot(Freq1(a-1,:),i,'.k'); hold on;
        
    end
    %     freqDIFF(j)=(mean(mean(ys(:,212:226))))-(mean(mean(ys(:,126:151))));
    freqDIFFfish2(j)=(mean(mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    perffish2(j)=sum((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5)/length((mean(ys(:,176:201),2)-mean(ys(:,126:151),2)));
    NRtrialsfish1=((mean(ys(:,176:201),2)-mean(ys(:,126:151),2))>=0.5);
    %     subplot 144; [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k'); xlim([-1 2.5])
    a=1; Interv=[];
    for i=1:size(sample4,1)
        [idx2,idx]=min(abs(sample3-sample4(i)));
        
        if sample3(idx)-sample4(i)<0
            idx=idx+1;
        end
        
        try
            Interv(a)=(sample3(idx)-sample4(i))/samplerate1;
        catch
            Interv(a)=0;
            %plot(Freq1(a-1,:),i,'.k'); hold on;
        end
        if Interv(a)<=0
            Interv(a)=0;
        end
        a=a+1;
    end
    %figure; hist(Interv,500)
    AUX=[]; AUX2=[]; AUX2=sample4(Interv>=0.010 & Interv<=0.020)/samplerate1; AUX3=[];
    for l=1:size(sample1_N,1)
        AUX3(l)=sum(AUX2>=sample1_N(l)/samplerate1 & AUX2<=(sample1_N(l)/samplerate1)+1);
    end
    AUX3(AUX3==0)=[];
    perfInterv(j)=size(AUX3,2)/size(sample1_N,1);
    
    
end

figure; [param,stat1]=sigm_fit(objdist,perfInterv,[nan nan nan -0.2]);
R2_Interv=1-(sum((stat1.ypred-perfInterv').^2)/sum((mean(perfInterv)-perfInterv).^2));


% figure; [param]=sigm_fit(objdist,freqDIFF,[NaN NaN 31 -0.2]);
% figure; [param,stat1]=sigm_fit(objdist,perf,[nan nan nan -0.2]);
% % R2_1=1-(sum((stat1.ypred-perf').^2)/sum((mean(perf)-perf).^2));
% figure; [param]=sigm_fit(objdist,freqDIFFfish2,[NaN NaN 31 -0.2]);
%
% figure; [param,stat2]=sigm_fit(objdist,perffish2,[nan nan nan -0.2]);
% R2_2=1-(sum((stat2.ypred-perffish2').^2)/sum((mean(perffish2)-perffish2).^2));

%% test 20-05-2021

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\Fede\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% second part, 2 fish

%objdist=[0 51 60.8 31.6 36.4 26.9 20 22.36 23.25 41.23 10 30];
objdist=[60.8 20 31.6 41.23 15 10 51 26.9];
vidnum=[15 17 18 19 20 21 22 25];

samplerate1=20000;
eventsamplerate=1000;
perfInterv=[];
 j=2;%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);  events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']); MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00003,'MINPEAKDISTANCE',50);
    if j>=3 && j<=4
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.005,'MINPEAKDISTANCE',50);
    else
        [value4,sample4]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.005,'MINPEAKDISTANCE',50);
    end
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
    
    for k=1:size(sample4,1)
        value3(sample3(:)==sample4(k))=[];
        sample3(sample3(:)==sample4(k))=[];
    end
    
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',4,'MINPEAKDISTANCE',1900);
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1; EODr1=zscore(EODr1);
    EODrate2=(diff(sample4)/samplerate1); EODr2=1./EODrate2;  EODr2(2:end+1)=EODr2; EODr2=zscore(EODr2);
    
    

    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+1)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+1))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])
    
    figure; subplot 121
    a=1; xs=[]; ys=[]; xt=[]; yt=[]; Freq1post=[]; Freq1=[];
    for i=1:size(sample1_N)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4(idx)>sample1_N(i)
            [idx3,idx4]=min(abs(sample3-sample4(idx+2)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+2))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        else
            [idx3,idx4]=min(abs(sample3-sample4(idx+3)));
            try
                Freq1post(a,:)=(sample3(idx4-40:idx4+40)-sample4(idx+3))./samplerate1;
                Freq1(a,:)=EODr1(idx4-40:idx4+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.999999999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end
        end
    end
    subplot 122; [hl, hp]=boundedline(xs', mean(ys),std(ys)/size(ys,1),'-k'); xlim([-1 2.5])