addpath('D:\KIT3'); addpath('C:\Users\fedu1\Google Drive\training')
clearvars; %close all;
myKsDir = uigetdir('F:\');
files2=dir([myKsDir, '\EODdata2*']);
a=1; b=1; xs=[]; ys=[]; Freq1post=[]; Freq2post=[]; Freq1=[]; samplerate1=30000;
%% obj as 0
tic
for j=1:length(files2)
    fileID = fopen([myKsDir,'\',files2(j).name]);
    A = fread(fileID,[4,Inf],'int16'); fclose(fileID);
    disp(['starting file ' num2str(j)])
    toc
    figure; plot(A(4,1:1000000)); 
    [~,yinput] = ginput(1); close
    [~,sample3]=findpeaks(A(4,:),'MINPEAKHEIGHT',yinput,'MINPEAKDISTANCE',200); %real fish + mimic
    [~,sample4]=findpeaks(A(3,:) ,'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',200); % just mimic
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime); hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
     
    for k=1:size(sample4,2) %remove mimic from fish + mimic
        sample3(sample3(:)>=sample4(k)-10 & sample3(:)<=sample4(k)+10)=[];
    end
    %% remove artefact
    [~,sample1_N]=findpeaks((zscore(A(1,:))) ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',60000); % find event, object-on
    for k=1:size(sample1_N,2) %remove mimic from fish + mimic
        sample3(sample3(:)>=sample1_N(k)-200 & sample3(:)<=sample1_N(k)+200)=[];
    end
    
    sample1_N=sample1_N-2350;
    
    for k=1:size(sample1_N,2) %remove mimic from fish + mimic
        sample3(sample3(:)>=sample1_N(k)-200 & sample3(:)<=sample1_N(k)+200)=[];
    end
    %%
   % figure; plot(1/samplerate1:1/samplerate1:size(A,2)/samplerate1,(A(1,:))); hold on; plot(sample1_N/samplerate1,1,'ok');    
    clear A    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); EODr1=zscore(EODr1);
    
    for i=1:size(sample1_N,2)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        try
            Freq1post(a,:)=(sample3(idx-120:idx+120)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-120:idx+120);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-60000/samplerate1 60000/samplerate1],0.9999999999);
            %plot(xs,ys(a,:),'-k'); hold on;
            a=a+1;
        end
    end    
    
    
    for i=1:size(sample1_N,2)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        try
            Freq2post(b,:)=(sample4(idx-120:idx+120)-sample1_N(i))./samplerate1;
            %plot(xs,ys(a,:),'-k'); hold on;
            b=b+1;
        end
    end
end
%%
figure; subplot(2,2,1); plot(xs(1,:),ys,'-k'); 
subplot(2,2,2);[hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');
subplot(2,2,3);
for i=1:size(Freq1post,1)
    plot(Freq1post(i,:),i,'.k'); xlim([-1 2]); hold on;
end
subplot(2,2,4);    
for i=1:size(Freq2post,1)
    plot(Freq2post(i,:),i,'.k'); xlim([-1 2]); hold on;
end    


%% fish as 0
tic
for j=1:length(files2)
    fileID = fopen([myKsDir,'\',files2(j).name]);
    A = fread(fileID,[4,Inf],'int16'); fclose(fileID);
    disp(['starting file ' num2str(j)])
    toc
    figure; plot(A(4,1:1000000)); 
    [~,yinput] = ginput(1); close
    [~,sample3]=findpeaks(A(4,:),'MINPEAKHEIGHT',yinput,'MINPEAKDISTANCE',200); %real fish + mimic
    [~,sample4]=findpeaks(A(3,:) ,'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',200); % just mimic
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime); hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
     
    for k=1:size(sample4,2) %remove mimic from fish + mimic
        sample3(sample3(:)>=sample4(k)-10 & sample3(:)<=sample4(k)+10)=[];
    end
    
    [~,sample1_N]=findpeaks((zscore(A(1,:))) ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',60000); % find event, object-on
    sample1_N=sample1_N-9000;
    figure; plot(1/samplerate1:1/samplerate1:size(A,2)/samplerate1,(A(1,:))); hold on; plot(sample1_N/samplerate1,1,'ok');    
    clear A    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); EODr1=zscore(EODr1);
    
    for i=1:size(sample1_N,2)
        [idx2,idx]=min(abs(sample3-sample1_N(i)));
        if sample3-sample1_N(i)>=0
            idx=idx;
        else
            idx=idx+1;
        end
        try
            Freq1post(a,:)=(sample3(idx-120:idx+120)-sample3(idx))./samplerate1;
            Freq1(a,:)=EODr1(idx-120:idx+120);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-60000/samplerate1 60000/samplerate1],0.9999999999);
            %plot(xs,ys(a,:),'-k'); hold on;
            a=a+1;
        end
    end    
    
    
    for i=1:size(sample1_N,2)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4-sample1_N(i)>=0
            idx=idx;
        else
            idx=idx+1;
        end      
        try
            Freq2post(b,:)=(sample4(idx-120:idx+120)-sample1_N(i))./samplerate1;
            %plot(xs,ys(a,:),'-k'); hold on;
            b=b+1;
        end
    end
end
%%
figure; subplot(2,2,1); plot(xs(1,:),ys,'-k'); 
subplot(2,2,2);[hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');
subplot(2,2,3);
for i=1:size(Freq1post,1)
    plot(Freq1post(i,:),i,'.k'); xlim([-1 2]); hold on;
end
subplot(2,2,4);    
for i=1:size(Freq2post,1)
    plot(Freq2post(i,:),i,'.k'); xlim([-1 2]); hold on;
end    
%% mimic as 0
tic
for j=1:length(files2)
    fileID = fopen([myKsDir,'\',files2(j).name]);
    A = fread(fileID,[4,Inf],'int16'); fclose(fileID);
    disp(['starting file ' num2str(j)])
    toc
    figure; plot(A(4,1:1000000)); 
    [~,yinput] = ginput(1); close
    [~,sample3]=findpeaks(A(4,:),'MINPEAKHEIGHT',yinput,'MINPEAKDISTANCE',200); %real fish + mimic
    [~,sample4]=findpeaks(A(3,:) ,'MINPEAKHEIGHT',0,'MINPEAKDISTANCE',200); % just mimic
    %figure; subplot(2,1,1); plot((1:1:size(EODtime,1))/samplerate1,EODtime); hold on; plot(sample3/samplerate1,value3,'or'); hold on; plot(sample4/samplerate1,value4,'ob');
     
    for k=1:size(sample4,2) %remove mimic from fish + mimic
        sample3(sample3(:)>=sample4(k)-10 & sample3(:)<=sample4(k)+10)=[];
    end
    
    [~,sample1_N]=findpeaks((zscore(A(1,:))) ,'MINPEAKHEIGHT',10,'MINPEAKDISTANCE',60000); % find event, object-on
    sample1_N=sample1_N-9000;
    %figure; plot(1/samplerate1:1/samplerate1:size(A,2)/samplerate1,zscore(A(1,:))); hold on; plot(sample1_N/samplerate1,1,'ok');    
    clear A    
    EODrate1=(diff(sample3)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); EODr1=zscore(EODr1);
    
    for i=1:size(sample1_N,2)
        [idx2,idx3]=min(abs(sample4-sample1_N(i)));
        if sample4-sample1_N(i)>=0
            idx3=idx3;
        else
            idx3=idx3+1;
        end
        [idx2,idx]=min(abs(sample3-sample4(idx3)));
        
        try
            Freq1post(a,:)=(sample3(idx-120:idx+120)-sample1_N(i))./samplerate1;
            Freq1(a,:)=EODr1(idx-120:idx+120);
            [xs(a,:), ys(a,:)]=FitVal_EI(Freq1post(a,:), Freq1(a,:), [-60000/samplerate1 60000/samplerate1],0.9999999999);
            %plot(xs,ys(a,:),'-k'); hold on;
            a=a+1;
        end
    end    
    
    
    for i=1:size(sample1_N,2)
        [idx2,idx]=min(abs(sample4-sample1_N(i)));
        if sample4-sample1_N(i)>=0
            idx=idx;
        else
            idx=idx+1;
        end      
        try
            Freq2post(b,:)=(sample4(idx-120:idx+120)-sample1_N(i))./samplerate1;
            %plot(xs,ys(a,:),'-k'); hold on;
            b=b+1;
        end
    end
end
%%
figure; subplot(2,2,1); plot(xs(1,:),ys,'-k'); 
subplot(2,2,2);[hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');
subplot(2,2,3);
for i=1:size(Freq1post,1)
    plot(Freq1post(i,:),i,'.k'); xlim([-1 2]); hold on;
end
subplot(2,2,4);    
for i=1:size(Freq2post,1)
    plot(Freq2post(i,:),i,'.k'); xlim([-1 2]); hold on;
end    
%figure; subplot(1,4,2);boxplot(EODr1,'PlotStyle','compact');
%subplot(1,4,1); [hl, hp]=boundedline(xs', mean(ys),std(ys),'-k');