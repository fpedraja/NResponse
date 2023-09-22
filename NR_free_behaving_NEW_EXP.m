addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\8Fish_new_exp_2\');
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

