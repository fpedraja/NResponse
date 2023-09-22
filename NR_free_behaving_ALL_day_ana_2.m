% FIRST PART
addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\6Fish_new_exp\');
files2=dir([myKsDir, '\EODdata2*']);
files3=dir([myKsDir, '\obj_num*']);
files4=dir([myKsDir, '\video*.avi']);
%%
a=1; b=1; c=1; xs=[]; ys=[]; Freq1post=[]; Freq2post=[]; Freq1=[]; samplerate1=30000;  FreqTotal_1=[]; FreqTotal_2=[]; Freqtotalpost_1=[]; Freqtotalpost_2=[]; FreqTotal_control=[]; Freqtotalpost__control=[];
Obj_idx_control=[]; Obj_idx_1=[]; Obj_idx_2=[]; M1total=[]; M2total=[]; Objtotal1=[]; Objtotal2=[];
 
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
            
                M1=[]; M2=[]; OB1=[]; OB2=[];
            if  size(sample1_N,2)==size(M,1)
                if abs(min(zscore(A(6,:))))<=abs(max(zscore(A(6,:))))
                    [~,sample1_FM1M2]=findpeaks((A(6,:)-min(A(6,:)))/ range(A(6,:)) ,'MINPEAKHEIGHT',0.55,'MINPEAKDISTANCE',100); % recording electrodes; fish + Mimic 1 + Mimic 2
                else
                    [~,sample1_FM1M2]=findpeaks((-A(6,:)-min(-A(6,:)))/ range(A(6,:)) ,'MINPEAKHEIGHT',0.55,'MINPEAKDISTANCE',100);
                end
                %% move to no mimic or mimic part
                if (size(sample1_M1,2)<2000 || isempty(sample1_M1)==1) && (size(sample1_M2,2)<2000 || isempty(sample1_M2)==1)  % no mimic / control
                    
                elseif (isempty(sample1_M2)==1 || size(sample1_M2,2)<500) && isempty(sample1_M1)==0     % mimic present at mimic 1
                    for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
                    end
                    
                    for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
                    end
                    
                    EODrate1=(diff(sample1_FM1M2)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); EODr1=zscore(EODr1);
                    
                    a=1;
                    for k=1:size(sample1_FM1M2,2) %remove mimic from fish + mimic
                        try
                            [idx2,idx]=min(abs(sample1_M1-sample1_FM1M2(k)));
                            if sample1_M1(idx)>=sample1_FM1M2(k)
                                M1(a,1)=(sample1_M1(idx-1)-sample1_FM1M2(k))/samplerate1;
                                M1(a,2)=(sample1_M1(idx)-sample1_FM1M2(k))/samplerate1;
                                M1(a,3)=(sample1_FM1M2(k)+b)/samplerate1;
                            else
                                M1(a,1)=(sample1_M1(idx)-sample1_FM1M2(k))/samplerate1;
                                M1(a,2)=(sample1_M1(idx+1)-sample1_FM1M2(k))/samplerate1;
                                M1(a,3)=(sample1_FM1M2(k)+b)/samplerate1;
                            end
                            a=a+1;
                        end
                    end
                    OB1(:,1)=(sample1_N+b)/samplerate1;
                    b=b+size(A,2)+1;
                    
                    
                elseif (isempty(sample1_M1)==1 || size(sample1_M1,2)<500) && isempty(sample1_M2)==0    % mimic present at mimic 2
                    for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
                    end
                    
                    for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
                    end

                    EODrate1=(diff(sample1_FM1M2)/samplerate1); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end); EODr1=zscore(EODr1);
                    
                    a=1;
                    for k=1:size(sample1_FM1M2,2) %remove mimic from fish + mimic
                        try
                            [idx2,idx]=min(abs(sample1_M2-sample1_FM1M2(k)));
                            if sample1_M2(idx)>=sample1_FM1M2(k)
                                M2(a,1)=(sample1_M2(idx-1)-sample1_FM1M2(k))/samplerate1;
                                M2(a,2)=(sample1_M2(idx)-sample1_FM1M2(k))/samplerate1;
                                M2(a,3)=(sample1_FM1M2(k)+c)/samplerate1;
                            else
                                M2(a,1)=(sample1_M2(idx)-sample1_FM1M2(k))/samplerate1;
                                M2(a,2)=(sample1_M2(idx+1)-sample1_FM1M2(k))/samplerate1;
                                M2(a,3)=(sample1_FM1M2(k)+c)/samplerate1;
                            end
                            a=a+1;
                        end
                    end
                    OB2(:,1)=(sample1_N+c)/samplerate1;
                    c=c+size(A,2)+1;
                end
            end
        end
    clear A    
    M1total=[M1total;M1];
    M2total=[M2total;M2];
    Objtotal1=[Objtotal1; OB1];
    Objtotal2=[Objtotal2; OB2];
    end
end
%%
figure; subplot(3,2,1:2); plot(M2total(:,3),M2total(:,1),'.b',M2total(:,3),M2total(:,2),'.r'); hold on; line([Objtotal2 Objtotal2],[-0.1 0.1]); hold on; line([Objtotal2+1 Objtotal2+1],[-0.1 0.1]); ylim([-0.05 0.05]);
hold on; line([0 max(M2total(:,3))],[-0.010 -0.010]); hold on; line([0 max(M2total(:,3))],[-0.015 -0.015]);
hold on; line([0 max(M2total(:,3))],[0.010 0.010]); hold on; line([0 max(M2total(:,3))],[0.015 0.015]);
subplot(3,2,3:4); hist3(M2total(:,[3 1]),'Nbins',[round(max(M2total(:,3))/4,-3) 25],'CDataMode','auto','FaceColor','interp');shading interp; ylim([-0.05 0]); view(2);
subplot(3,2,5); hist(M2total(:,1),500); xlim([-0.05 0]); subplot(3,2,6); hist(M2total(:,2),500); xlim([0 0.05]);

figure; subplot(3,2,1:2); plot(M1total(:,3),M1total(:,1),'.b',M1total(:,3),M1total(:,2),'.r'); hold on; line([Objtotal1 Objtotal1],[-0.1 0.1]); hold on; line([Objtotal1+1 Objtotal1+1],[-0.1 0.1]); ylim([-0.05 0.05]);
hold on; line([0 max(M1total(:,3))],[-0.010 -0.010]); hold on; line([0 max(M1total(:,3))],[-0.015 -0.015]);
hold on; line([0 max(M1total(:,3))],[0.010 0.010]); hold on; line([0 max(M1total(:,3))],[0.015 0.015]);
subplot(3,2,3:4); hist3(M1total(:,[3 1]),'Nbins',[round(max(M1total(:,3))/4,-3) 25],'CDataMode','auto','FaceColor','interp');shading interp; ylim([-0.05 0])
subplot(3,2,5); hist(M1total(:,1),500); xlim([-0.05 0]); subplot(3,2,6); hist(M1total(:,2),500); xlim([0 0.05]);
