%% 04.06.21 1 second

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% first part, 1 second 
objdist=[20 31.6 60.8 0 10 51 41.23];
vidnum=[15 16 17 18 19 20 21];



samplerate1=20000;
eventsamplerate=1000;

for j=1%1:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.001,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',300);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',1900);
    
    figure;
    subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    subplot(3,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(round(size(sample2,1)/140),10); Freqfish_E=zeros(round(size(sample1,1)/120),10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 7 8 9 10]
            try
                Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+19)));
            catch
                Freqfish_M(t,k)=nan;
            end
            a=a+20;
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+20;
        end
    end
    
    Fish_median=nanmedian(Freqfish_M); figure; plot(Fish_median,'-ok');
    figure; boxplot(Freqfish_M,'PlotStyle','compact')
end

%% 07.06.21 1 second

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% first part, 1 second
objdist=[31.6 20 60.8 0 10 51 41.23];
vidnum=[1 2 3 4 5 6 7];

samplerate1=20000;
eventsamplerate=1000;

for j=2%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.001,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',300);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',1900);
    
    figure;
    subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    subplot(3,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(round(size(sample2,1)/140),10); Freqfish_E=zeros(round(size(sample1,1)/120),10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 7 8 9 10]
            try
                Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+19)));
            catch
                Freqfish_M(t,k)=nan;
            end
            a=a+20;
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+20;
        end
    end
    
    Fish_median=nanmedian(Freqfish_M); figure; plot(Fish_median,'-ok');
    figure; boxplot(Freqfish_M,'PlotStyle','compact')
end

%% 07.06.21 5 second

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[20 60.8 0 10 51 41.23 31.6 ];
vidnum=[9 10 11 12 13 14 15];

samplerate1=20000;
eventsamplerate=1000;

for j=2%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.001,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',300);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',1900);
    
    figure;
    subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    subplot(3,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate); sample2_N=sample2*(samplerate1/secondfish);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,10); Freqfish_E=zeros(1,10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 7 8 9 10]
            try
                Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+24)));
            catch
                Freqfish_M(t,k)=nan;
            end
            a=a+25;
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+24)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+25;
        end
    end
    figure; plot(Freqfish_M,'-ok');
    
    
    % change in freq
   CHAfish_M=zeros(25,10); a=1;
   c=1;
   for k=[1 2 3 7 8 9 10]
       for t=1:25
           try
               AUX1=[]; AUX2=[];
               AUX1=mean(EODr1(Fish_time>=Mfish_time(a)-c & Fish_time<Mfish_time(a)));
               AUX2=mean(EODr1(Fish_time>Mfish_time(a) & Fish_time<=Mfish_time(a)+c));
               CHAfish_M(t,k)=AUX2-AUX1;
           catch
               CHAfish_M(t,k)=nan;
           end
           a=a+1;
       end
   end
   
   a=1;
   for k=[1 2 3 4 5 6]
       for t=1:25
           if k>=4
               try
                   AUX1=[]; AUX2=[];
                   AUX1=mean(EODr1(Fish_time>=event_time(a)-c & Fish_time<event_time(a)));
                   AUX2=mean(EODr1(Fish_time>event_time(a) & Fish_time<=event_time(a)+c));
                   CHAfish_M(t,k)=AUX2-AUX1;
               catch
                   CHAfish_M(t,k)=nan;
               end
           end
           a=a+1;
       end
   end
   
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
    
        a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; 

        for i=100:125
            [idx2,idx]=min(abs(sample3-sample1_N(i)));
            try
                Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end

  a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
    figure; 

        for i=101:125
            [idx2,idx]=min(abs(sample3-sample2_N(i)));
            try
                Freq1post(a,:)=(sample3(idx-40:idx+40)-sample2_N(i))./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end

    
    
    
end

AUX3=sum(CHAfish_M>0.4)/size(CHAfish_M,1);
%% 08.06.21 1 second

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% first part, 1 second
objdist=[31.6 20 60.8 0 10 51 41.23];
vidnum=[1 2 3 4 5 6 7];

samplerate1=20000;
eventsamplerate=1000;

for j=3%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.001,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',300);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',1900);
    
    figure;
    subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    subplot(3,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(round(size(sample2,1)/140),10); Freqfish_E=zeros(round(size(sample1,1)/120),10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 7 8 9 10]
            try
                Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+19)));
            catch
                Freqfish_M(t,k)=nan;
            end
            a=a+20;
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+20;
        end
    end
    
    Fish_median=nanmedian(Freqfish_M); figure; plot(Fish_median,'-ok');
    figure; boxplot(Freqfish_M,'PlotStyle','compact')
end

%% 08.06.21 5 second

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[0 60.8 41.23  31.6 10 51 20]; %(+10)
vidnum=[8 9 10 11 12 13 14];

samplerate1=20000;
eventsamplerate=1000;

for j=2%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0015,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',300);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',1900);
    
    figure;
    subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    subplot(3,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,10); Freqfish_E=zeros(1,10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 7 8 9 10]
            try
                Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+24)));
            catch
                Freqfish_M(t,k)=nan;
            end
            a=a+25;
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+24)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+25;
        end
    end
    figure; plot(Freqfish_M,'-ok');
    
    
    % change in freq
   CHAfish_M=zeros(25,10); a=1;
   for k=[1 2 3 7 8 9 10]
       for t=1:25
           try
               AUX1=[]; AUX2=[];
               AUX1=mean(EODr1(Fish_time>=Mfish_time(a)-1 & Fish_time<=Mfish_time(a)));
               AUX2=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a)+1));
               CHAfish_M(t,k)=AUX2-AUX1;
           catch
               CHAfish_M(t,k)=nan;
           end
           a=a+1;
       end
   end
   
   a=1;
   for k=[1 2 3 4 5 6]
       for t=1:25
           if k>=4
                try
                   AUX1=[]; AUX2=[];
                   AUX1=mean(EODr1(Fish_time>=event_time(a)-1 & Fish_time<=event_time(a)));
                   AUX2=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a)+1));
                   CHAfish_M(t,k)=AUX2-AUX1;
                catch
                    CHAfish_M(t,k)=nan;
                end
           end
           a=a+1;
       end
   end
   
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
end

AUX3=sum(CHAfish_M>0.4)/size(CHAfish_M,1);
%% 10.06.21 1 second

clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

% first part, 1 second
objdist=[31.6 20 60.8 0 10 51 41.23];
vidnum=[1 2 3 4 5 6 7];

samplerate1=20000;
eventsamplerate=1000;

for j=2%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.001,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',300);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',1900);
    
    figure;
    subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    subplot(3,1,2);  plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
    subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(round(size(sample2,1)/140),10); Freqfish_E=zeros(round(size(sample1,1)/120),10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 7 8 9 10]
            try
                Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+19)));
            catch
                Freqfish_M(t,k)=nan;
            end
            a=a+20;
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+20;
        end
    end
    
    Fish_median=nanmedian(Freqfish_M); figure; plot(Fish_median,'-ok');
    figure; boxplot(Freqfish_M,'PlotStyle','compact')
end

%% 10.06.21 5 second
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[31.6 20 60 30 20]; %
vidnum=[1 2 3 4 5];

samplerate1=20000;
eventsamplerate=1000;

for j=3%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0006,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',300);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',1900);
%     
     figure;
     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
    subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,10); Freqfish_E=zeros(1,10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 7 8 9 10]
            try
                Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+24)));
            catch
                Freqfish_M(t,k)=nan;
            end
            a=a+25;
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+24)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+25;
        end
    end
%     figure; plot(Freqfish_M,'-ok');
    
    c=1;
    % change in freq
   CHAfish_M=zeros(25,10); a=1;
   for k=[1 2 3 7 8 9 10]
       for t=1:25
           try
               AUX1=[]; AUX2=[];
               AUX1=mean(EODr1(Fish_time>=Mfish_time(a)-c & Fish_time<=Mfish_time(a)));
               AUX2=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a)+c));
               CHAfish_M(t,k)=AUX2-AUX1;
           catch
               CHAfish_M(t,k)=nan;
           end
           a=a+1;
       end
   end
   
   a=1;
   for k=[1 2 3 4 5 6]
       for t=1:25
           if k>=4
                try
                   AUX1=[]; AUX2=[];
                   AUX1=mean(EODr1(Fish_time>=event_time(a)-c & Fish_time<=event_time(a)));
                   AUX2=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a)+c));
                   CHAfish_M(t,k)=AUX2-AUX1;
                catch
                    CHAfish_M(t,k)=nan;
                end
           end
           a=a+1;
       end
   end
   
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
end
AUX3=sum(CHAfish_M>0.4)/size(CHAfish_M,1);


%% 11.06.21 5 second
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[60 0 20 31.6 20 30]; %
vidnum=[1 2 3 4 5 6];

samplerate1=20000;
eventsamplerate=1000;

for j=2%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0015,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',300);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',1900);
%     
%     figure;
%     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
%     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
%     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,10); Freqfish_E=zeros(1,10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 7 8 9 10]
            try
                Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+24)));
            catch
                Freqfish_M(t,k)=nan;
            end
            a=a+25;
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+24)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+25;
        end
    end
%     figure; plot(Freqfish_M,'-ok');
    
    c=2;
    % change in freq
   CHAfish_M=zeros(25,10); a=1;
   for k=[1 2 3 7 8 9 10]
       for t=1:25
           try
               AUX1=[]; AUX2=[];
               AUX1=mean(EODr1(Fish_time>=Mfish_time(a)-c & Fish_time<=Mfish_time(a)));
               AUX2=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a)+c));
               CHAfish_M(t,k)=AUX2-AUX1;
           catch
               CHAfish_M(t,k)=nan;
           end
           a=a+1;
       end
   end
   
   a=1;
   for k=[1 2 3 4 5 6]
       for t=1:25
           if k>=4
                try
                   AUX1=[]; AUX2=[];
                   AUX1=mean(EODr1(Fish_time>=event_time(a)-c & Fish_time<=event_time(a)));
                   AUX2=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a)+c));
                   CHAfish_M(t,k)=AUX2-AUX1;
                catch
                    CHAfish_M(t,k)=nan;
                end
           end
           a=a+1;
       end
   end
   
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
end

AUX3=sum(CHAfish_M>0.4)/size(CHAfish_M,1);


%% 18.06.21 5 second
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[10 20 30 60]; %
vidnum=[1 2 3 4 ];

samplerate1=20000;
eventsamplerate=1000;

for j=2%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0015,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',300);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',1900);
%     
%     figure;
%     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
%     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
%     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,10); Freqfish_E=zeros(1,10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 7 8 9 10]
            try
                Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+24)));
            catch
                Freqfish_M(t,k)=nan;
            end
            a=a+25;
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+24)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+25;
        end
    end
%     figure; plot(Freqfish_M,'-ok');
    
    c=3;
    % change in freq
   CHAfish_M=zeros(25,10); a=1;
   for k=[1 2 3 7 8 9 10]
       for t=1:25
           try
               AUX1=[]; AUX2=[];
               AUX1=mean(EODr1(Fish_time>=Mfish_time(a)-c & Fish_time<=Mfish_time(a)));
               AUX2=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a)+c));
               CHAfish_M(t,k)=AUX2-AUX1;
           catch
               CHAfish_M(t,k)=nan;
           end
           a=a+1;
       end
   end
   
   a=1;
   for k=[1 2 3 4 5 6]
       for t=1:25
           if k>=4
                try
                   AUX1=[]; AUX2=[];
                   AUX1=mean(EODr1(Fish_time>=event_time(a)-c & Fish_time<=event_time(a)));
                   AUX2=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a)+c));
                   CHAfish_M(t,k)=AUX2-AUX1;
                catch
                    CHAfish_M(t,k)=nan;
                end
           end
           a=a+1;
       end
   end
   
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
end

AUX3=sum(CHAfish_M>0.4)/size(CHAfish_M,1)

%% 18.06.21 20 seconds
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[60 30 20 10 20]; %
vidnum=[5 6 7 8 9];

samplerate1=20000;
eventsamplerate=1000;

for j=1%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0015,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',2000);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',300);
%     
     figure;
     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,10); Freqfish_E=zeros(1,10);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6 7 8 9 10]
            if k>=4 && k<=6
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+377)));
                catch
                    Freqfish_M(t,k)=nan;
                end
                a=a+378;
            else
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a+440)));
                catch
                    Freqfish_M(t,k)=nan;
                end
                a=a+441;
            end
        end
    end
    
    %      a=1;
    %    for t=1:size(Freqfish_E,1)
    %        for k=[1 2 3 4 5 6]
    %         Freqfish_E(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+19)));
    %         a=a+20;
    %        end
    %    end
    
    
    a=1;
    for t=1:size(Freqfish_M,1)
        for k=[1 2 3 4 5 6]
            if k>=4
                try
                    Freqfish_M(t,k)=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a+20)));
                catch
                    Freqfish_M(t,k)=nan;
                end
            end
            a=a+21;
        end
    end
%     figure; plot(Freqfish_M,'-ok');
%%    
    c=1;
    % change in freq
   CHAfish_M=zeros(21,10); 
   a=1;
   for k=[1 2 3 4 5 6 7 8 9 10]
       if k>=4 && k<=6
           for t=1:21
               try
                   AUX1=[]; AUX2=[];
                   AUX1=mean(EODr1(Fish_time>=Mfish_time(a)-c & Fish_time<=Mfish_time(a)));
                   AUX2=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a)+c));
                   CHAfish_M(t,k)=AUX2-AUX1;
               catch
                   CHAfish_M(t,k)=nan;
               end
               a=a+18;
           end
       else
           for t=1:21
               try
                   AUX1=[]; AUX2=[];
                   AUX1=mean(EODr1(Fish_time>=Mfish_time(a)-c & Fish_time<=Mfish_time(a)));
                   AUX2=mean(EODr1(Fish_time>=Mfish_time(a) & Fish_time<=Mfish_time(a)+c));
                   CHAfish_M(t,k)=AUX2-AUX1;
               catch
                   CHAfish_M(t,k)=nan;
               end
               a=a+21;
           end
       end
   end
   
%    a=1;
%    for k=[1 2 3 4 5 6]
%        for t=1:25
%            if k>=4
%                 try
%                    AUX1=[]; AUX2=[];
%                    AUX1=mean(EODr1(Fish_time>=event_time(a)-c & Fish_time<=event_time(a)));
%                    AUX2=mean(EODr1(Fish_time>=event_time(a) & Fish_time<=event_time(a)+c));
%                    CHAfish_M(t,k)=AUX2-AUX1;
%                 catch
%                     CHAfish_M(t,k)=nan;
%                 end
%            end
%            a=a+1;
%        end
%    end
   
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
    
     a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
     D=[1    22    43    64    85   106 ];
    
for r=1:6
    a=1; figure; 
        for i=D(r):D(r)+21
            [idx2,idx]=min(abs(sample3-sample1_N(i)));
            try
                Freq1post(a,:)=(sample3(idx-40:idx+40)-sample1_N(i))./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end 
end
end

AUX3=sum(CHAfish_M>0.8)/size(CHAfish_M,1)

%% 22.06.21 20 seconds
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[30 20 60 30]; %
vidnum=[1 2 3 4];

samplerate1=20000;
eventsamplerate=1000;

for j=3%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0025,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',2000);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',300);
%     
     figure;
     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,11); Freqfish_E=zeros(1,11);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
%%
 load('C:\Users\fedu1\Desktop\positions.mat')  
   
 c=3;
 % change in freq
 CHAfish_M=zeros(21,11);
 for i=1:size(positions,1)
     for j=1:size(positions,2)        
         try
             AUX1=[]; AUX2=[];
             AUX1=mean(EODr1(Fish_time>=Mfish_time(positions(i,j))-c & Fish_time<=Mfish_time(positions(i,j))));
             AUX2=mean(EODr1(Fish_time>=Mfish_time(positions(i,j)) & Fish_time<=Mfish_time(positions(i,j))+c));
             CHAfish_M(i,j)=AUX2-AUX1;
         catch
             CHAfish_M(i,j)=nan;
         end
     end
 end
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
    
%      a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
%      D=[1    22    43    64    85   106 ];
%     

    a=1; figure; 
        for i=1:21
            [idx2,idx]=min(abs(Fish_time-Mfish_time(positions(i,11))));
            try
                Freq1post(a,:)=(Fish_time(idx-40:idx+40)-Mfish_time(positions(i,11)));%./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end 

end

AUX3=sum(CHAfish_M>0.6)/size(CHAfish_M,1)

%% 24.06.21 20 seconds
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[30 20 60 30]; %
vidnum=[1 2 3 4 5 6];

samplerate1=20000;
eventsamplerate=1000;

for j=4%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0025,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',2000);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',300);
%     
     figure;
     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,11); Freqfish_E=zeros(1,11);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
%%
 load('C:\Users\fedu1\Desktop\positions.mat')  
   
 c=2.5;
 % change in freq
 CHAfish_M=zeros(21,11);
 for i=1:size(positions,1)
     for j=1:size(positions,2)        
         try
             AUX1=[]; AUX2=[];
             AUX1=mean(EODr1(Fish_time>=Mfish_time(positions(i,j))-c & Fish_time<=Mfish_time(positions(i,j))));
             AUX2=mean(EODr1(Fish_time>=Mfish_time(positions(i,j)) & Fish_time<=Mfish_time(positions(i,j))+c));
             CHAfish_M(i,j)=AUX2-AUX1;
         catch
             CHAfish_M(i,j)=nan;
         end
     end
 end
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
    
%      a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
%      D=[1    22    43    64    85   106 ];
%     

    a=1; figure; 
        for i=1:21
            [idx2,idx]=min(abs(Fish_time-Mfish_time(positions(i,11))));
            try
                Freq1post(a,:)=(Fish_time(idx-40:idx+40)-Mfish_time(positions(i,11)));%./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end 

end

AUX3=sum(CHAfish_M>0.5)/size(CHAfish_M,1)

%% 25.06.21 20 seconds
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[30 20 60 30]; %
vidnum=[1 2 3 4 5 6];

samplerate1=20000;
eventsamplerate=1000;

for j=4%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0025,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',2000);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',300);
%     
     figure;
     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,11); Freqfish_E=zeros(1,11);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
%%
 load('C:\Users\fedu1\Desktop\positions.mat')  
   
 c=2.5;
 % change in freq
 CHAfish_M=zeros(21,11);
 for i=1:size(positions,1)
     for j=1:size(positions,2)        
         try
             AUX1=[]; AUX2=[];
             AUX1=mean(EODr1(Fish_time>=Mfish_time(positions(i,j))-c & Fish_time<=Mfish_time(positions(i,j))));
             AUX2=mean(EODr1(Fish_time>=Mfish_time(positions(i,j)) & Fish_time<=Mfish_time(positions(i,j))+c));
             CHAfish_M(i,j)=AUX2-AUX1;
         catch
             CHAfish_M(i,j)=nan;
         end
     end
 end
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
    
%      a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
%      D=[1    22    43    64    85   106 ];
%     

    a=1; figure; 
        for i=1:21
            [idx2,idx]=min(abs(Fish_time-Mfish_time(positions(i,11))));
            try
                Freq1post(a,:)=(Fish_time(idx-40:idx+40)-Mfish_time(positions(i,11)));%./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end 

end

AUX3=sum(CHAfish_M>0.5)/size(CHAfish_M,1)


%% 29.06.21 20 seconds
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[30 20 60 30]; %
vidnum=[1 2 3 4 5];

samplerate1=20000;
eventsamplerate=1000;

for j=4%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00028,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',2000);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',300);
%     
     figure;
     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,11); Freqfish_E=zeros(1,11);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
%%
 load('C:\Users\fedu1\Desktop\positions.mat')  
   
 c=2.5;
 % change in freq
 CHAfish_M=zeros(21,11);
 for i=1:size(positions,1)
     for j=1:size(positions,2)        
         try
             AUX1=[]; AUX2=[];
             AUX1=mean(EODr1(Fish_time>=Mfish_time(positions(i,j))-c & Fish_time<=Mfish_time(positions(i,j))));
             AUX2=mean(EODr1(Fish_time>=Mfish_time(positions(i,j)) & Fish_time<=Mfish_time(positions(i,j))+c));
             CHAfish_M(i,j)=AUX2-AUX1;
         catch
             CHAfish_M(i,j)=nan;
         end
     end
 end
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
    
%      a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
%      D=[1    22    43    64    85   106 ];
%     

    a=1; figure; 
        for i=1:21
            [idx2,idx]=min(abs(Fish_time-Mfish_time(positions(i,11))));
            try
                Freq1post(a,:)=(Fish_time(idx-40:idx+40)-Mfish_time(positions(i,11)));%./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end 

end

AUX3=sum(CHAfish_M>0.6)/size(CHAfish_M,1)

%% 30.06.21 20 seconds
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[30 20 60 30]; %
vidnum=[1 2 3 4 5];

samplerate1=20000;
eventsamplerate=1000;

for j=3%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0005,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',2000);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',300);
%     
     figure;
     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,11); Freqfish_E=zeros(1,11);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
%%
 load('C:\Users\fedu1\Desktop\positions.mat')  
   
 c=2.5;
 % change in freq
 CHAfish_M=zeros(21,11);
 for i=1:size(positions,1)
     for j=1:size(positions,2)        
         try
             AUX1=[]; AUX2=[];
             AUX1=mean(EODr1(Fish_time>=Mfish_time(positions(i,j))-c & Fish_time<=Mfish_time(positions(i,j))));
             AUX2=mean(EODr1(Fish_time>=Mfish_time(positions(i,j)) & Fish_time<=Mfish_time(positions(i,j))+c));
             CHAfish_M(i,j)=AUX2-AUX1;
         catch
             CHAfish_M(i,j)=nan;
         end
     end
 end
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
    
%      a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
%      D=[1    22    43    64    85   106 ];
%     

    a=1; figure; 
        for i=1:21
            [idx2,idx]=min(abs(Fish_time-Mfish_time(positions(i,11))));
            try
                Freq1post(a,:)=(Fish_time(idx-40:idx+40)-Mfish_time(positions(i,11)));%./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end 

end

AUX3=sum(CHAfish_M>0.6)/size(CHAfish_M,1)

%% 01.07.21 20 seconds
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[30 20 60 30]; %
vidnum=[1 2 3 4 5];

samplerate1=20000;
eventsamplerate=1000;

for j=2%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.00079,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',2000);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',300);
%     
     figure;
     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,11); Freqfish_E=zeros(1,11);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
%%
 load('C:\Users\fedu1\Desktop\positions.mat')  
   
 c=2.5;
 % change in freq
 CHAfish_M=zeros(21,11);
 for i=1:size(positions,1)
     for j=1:size(positions,2)        
         try
             AUX1=[]; AUX2=[];
             AUX1=mean(EODr1(Fish_time>=Mfish_time(positions(i,j))-c & Fish_time<=Mfish_time(positions(i,j))));
             AUX2=mean(EODr1(Fish_time>=Mfish_time(positions(i,j)) & Fish_time<=Mfish_time(positions(i,j))+c));
             CHAfish_M(i,j)=AUX2-AUX1;
         catch
             CHAfish_M(i,j)=nan;
         end
     end
 end
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
    
%      a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
%      D=[1    22    43    64    85   106 ];
%     

    a=1; figure; 
        for i=1:21
            [idx2,idx]=min(abs(Fish_time-Mfish_time(positions(i,11))));
            try
                Freq1post(a,:)=(Fish_time(idx-40:idx+40)-Mfish_time(positions(i,11)));%./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end 

end

AUX3=sum(CHAfish_M>0.6)/size(CHAfish_M,1)

%% 02.07.21 20 seconds
clearvars; cd('D:\KIT'); addpath('D:\KIT3'); addpath('D:\KIT2');
direfinal=uigetdir('Z:\locker\2FISH\'); %3D data
files2=dir([direfinal, '\*.mat']);

objdist=[30 20 60 30]; %
vidnum=[1 2 3 4 5];

samplerate1=20000;
eventsamplerate=1000;

for j=2%:length(vidnum)
    load([direfinal, '\', files2(vidnum(j)).name])
    EODtime=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch2.values']);
    events=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch4.values']);
    Mfish=eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.values']);
    secondfish=1/eval(['V',files2(vidnum(j)).name(1:end-4),'_Ch5.interval']);
    MINVAL1=[]; MINVAL2=[];
    
    [value3,sample3]=findpeaks(EODtime ,'MINPEAKHEIGHT',0.0009,'MINPEAKDISTANCE',90);
    [value1,sample1]=findpeaks(events ,'MINPEAKHEIGHT',3,'MINPEAKDISTANCE',2000);
    [value2,sample2]=findpeaks(Mfish ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',300);
%     
     figure;
     subplot(3,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime);  hold on; plot(sample3/samplerate1,value3,'or');
     subplot(3,1,2); plot((1:1:size(events,1))/eventsamplerate,events);  hold on; plot(sample1/eventsamplerate,value1,'ok');
     subplot(3,1,3); plot((1:1:size(Mfish,1))/secondfish,Mfish);  hold on; plot(sample2/secondfish,value2,'or');
    
    clear EODtime events Mfish
    
    sample1_N=sample1*(samplerate1/eventsamplerate);
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1=zscore(EODr1);
    
    Freqfish_M=zeros(1,11); Freqfish_E=zeros(1,11);
    Fish_time=(sample3)/samplerate1; Mfish_time=sample2/secondfish; event_time=sample1/eventsamplerate;
    
%%
 load('C:\Users\fedu1\Desktop\positions.mat')  
   
 c=2.5;
 % change in freq
 CHAfish_M=zeros(21,11);
 for i=1:size(positions,1)
     for j=1:size(positions,2)        
         try
             AUX1=[]; AUX2=[];
             AUX1=mean(EODr1(Fish_time>=Mfish_time(positions(i,j))-c & Fish_time<=Mfish_time(positions(i,j))));
             AUX2=mean(EODr1(Fish_time>=Mfish_time(positions(i,j)) & Fish_time<=Mfish_time(positions(i,j))+c));
             CHAfish_M(i,j)=AUX2-AUX1;
         catch
             CHAfish_M(i,j)=nan;
         end
     end
 end
   
    figure; subplot 311; boxplot(CHAfish_M,'PlotStyle','compact')
    subplot 312; plot(nanmean(CHAfish_M),'-ok')
    subplot 313; plot(nanmedian(CHAfish_M),'-ok')
    
%      a=1; xs=[]; ys=[]; Freq1post=[]; Freq1=[];
%      D=[1    22    43    64    85   106 ];
%     

    a=1; figure; 
        for i=1:21
            [idx2,idx]=min(abs(Fish_time-Mfish_time(positions(i,11))));
            try
                Freq1post(a,:)=(Fish_time(idx-40:idx+40)-Mfish_time(positions(i,11)));%./samplerate1;
                Freq1(a,:)=EODr1(idx-40:idx+40);
                [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-70000/samplerate1 70000/samplerate1],0.9999);
                plot(Freq1post(a,:),a,'.k'); xlim([-2 3.5]); hold on;
                a=a+1;
            end            
            %         plot(Freq1(a-1,:),i,'.k'); hold on;
        end 

end

AUX3=sum(CHAfish_M>0.6)/size(CHAfish_M,1)