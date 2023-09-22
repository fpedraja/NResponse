% FIRST PART
addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\1Fish_new_sil\20220927\');
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
            [~,sample1_N]=findpeaks((A(4,6e05:end)-min(A(4,6e05:end)))/ range(A(4,6e05:end)),'MINPEAKHEIGHT',0.8,'MINPEAKDISTANCE',90000); % find event, object-on
            
            
            if size(sample1_N,2)<size(M,1) %check if events in csv file corresponds to events detected
                [~,sample1_N]=findpeaks((A(4,6e05:end)-min(A(4,6e05:end)))/ range(A(4,6e05:end)) ,'MINPEAKHEIGHT',0.7,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            elseif size(sample1_N,2)>size(M,1) %check if events in csv file corresponds to events detected
                [~,sample1_N]=findpeaks((A(4,6e05:end)-min(A(4,6e05:end)))/ range(A(4,6e05:end)) ,'MINPEAKHEIGHT',0.9,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            end
            
            if size(sample1_N,2)<size(M,1) %check if events in csv file corresponds to events detected
                [~,sample1_N]=findpeaks((A(2,6e05:end)-min(A(2,6e05:end)))/ range(A(2,6e05:end)) ,'MINPEAKHEIGHT',0.7,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            elseif size(sample1_N,2)>size(M,1) %check if events in csv file corresponds to events detected
                [~,sample1_N]=findpeaks((A(2,6e05:end)-min(A(2,6e05:end)))/ range(A(2,6e05:end)) ,'MINPEAKHEIGHT',0.9,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            end
            
            if size(sample1_N,2)<size(M,1) %check if events in csv file corresponds to events detected
                [~,sample1_N]=findpeaks((A(1,6e05:end)-min(A(1,6e05:end)))/ range(A(1,6e05:end)) ,'MINPEAKHEIGHT',0.7,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            elseif size(sample1_N,2)>size(M,1) %check if events in csv file corresponds to events detected
                [~,sample1_N]=findpeaks((A(1,6e05:end)-min(A(1,6e05:end)))/ range(A(1,6e05:end)) ,'MINPEAKHEIGHT',0.9,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            end
            
            sample1_N=sample1_N+6e05;
            
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
                %% move to no mimic or mimic part
                if (size(sample1_M1,2)<2000 || isempty(sample1_M1)==1) && (size(sample1_M2,2)<2000 || isempty(sample1_M2)==1)  % no mimic / control
                    clear A
                    
                    
                    a=1; Freq1=[]; Freq1post=[]; AUX3=[]; AUX4=[];
                    vid = VideoReader([myKsDir,'\',files4(i).name]);
                    outputVideo = VideoWriter([myKsDir,'\CUT_',files4(i).name]);
                    outputVideo.FrameRate = vid.FrameRate;
                    open(outputVideo);
                    
                    for j=1:size(sample1_N,2)
                        
                       
                            for t=(round(sample1_N(j)/samplerate1*vid.FrameRate))-450:(round(sample1_N(j)/samplerate1*vid.FrameRate))+900
                                img = read(vid, t);  writeVideo(outputVideo,img);
                            end
                            AUX3(a,1)=M(j);
                            AUX4(a,1)=i;
                            a=a+1;
                        
                    end
                    Obj_idx_control=[Obj_idx_control;[AUX3 AUX4]];
                    close(outputVideo);
                    
                    
                    
                elseif (isempty(sample1_M2)==1 || size(sample1_M2,2)<500) && isempty(sample1_M1)==0     % mimic present at mimic 1
                    
                    clear A
                    a=1; Freq1=[]; Freq1post_1=[]; AUX3=[]; AUX4=[];
                    vid = VideoReader([myKsDir,'\',files4(i).name]);
                    outputVideo = VideoWriter([myKsDir,'\CUT_',files4(i).name]);
                    outputVideo.FrameRate = vid.FrameRate;
                    open(outputVideo);
                    
                    for j=1:size(sample1_N,2)
                        
                        
                            for t=(round(sample1_N(j)/samplerate1*vid.FrameRate))-450:(round(sample1_N(j)/samplerate1*vid.FrameRate))+900
                                img = read(vid, t);  writeVideo(outputVideo,img);
                            end
                            AUX3(a,1)=M(j);
                            AUX4(a,1)=i;
                            a=a+1;
                        
                    end
                    Obj_idx_1=[Obj_idx_1;[AUX3 AUX4]];
                    close(outputVideo);
                    
                    
                    
                    
                elseif (isempty(sample1_M1)==1 || size(sample1_M1,2)<500) && isempty(sample1_M2)==0    % mimic present at mimic 2
                    
                    
                    clear A
                    a=1; Freq1=[]; Freq1post_2=[]; AUX3=[]; AUX4=[];
                    vid = VideoReader([myKsDir,'\',files4(i).name]);
                    outputVideo = VideoWriter([myKsDir,'\CUT_',files4(i).name]);
                    outputVideo.FrameRate = vid.FrameRate;
                    open(outputVideo);
                    
                    for j=1:size(sample1_N,2)
                          
                            for t=(round(sample1_N(j)/samplerate1*vid.FrameRate))-450:(round(sample1_N(j)/samplerate1*vid.FrameRate))+900
                                img = read(vid, t);  writeVideo(outputVideo,img);
                            end
                            AUX3(a,1)=M(j);
                            AUX4(a,1)=i;
                            a=a+1;
                       
                    end
                    Obj_idx_2=[Obj_idx_2;[AUX3 AUX4]];
                    close(outputVideo);
                    
                    
                    
                end
            end
        end
    end
end

save([myKsDir(1:end-9),'\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_EOD_data.mat'], 'Obj_idx_1','Obj_idx_2', 'Obj_idx_control')

