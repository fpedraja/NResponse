% FIRST PART
addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\5Fish_new_exp\');
files2=dir([myKsDir, '\EODdata2*']);
files3=dir([myKsDir, '\obj_num*']);
files4=dir([myKsDir, '\video*.avi']);

a=1; b=1; c=1; xs=[]; ys=[]; Freq1post=[]; Freq2post=[]; Freq1=[]; samplerate1=30000;  FreqTotal_1=[]; FreqTotal_2=[]; Freqtotalpost_1=[]; Freqtotalpost_2=[]; FreqTotal_control=[]; Freqtotalpost__control=[];
Obj_idx_control=[]; Obj_idx_1=[]; Obj_idx_2=[]; M1=[]; M2=[];
%%
for i=1:size(files3,1)
    try %#ok<TRYNC> % check that the obj index has actual values (night vs day cycle)
        M = csvread([myKsDir,'\',files3(i).name]);
        if isempty(M)==1
        else
            fileID = fopen([myKsDir,'\',files2(i).name]);
            A = fread(fileID,[8,Inf],'int16'); fclose(fileID);
            %B = smooth(A(4,:),0.1,'moving');
            [~,sample1_M1]=findpeaks(((A(1,:))) ,'MINPEAKHEIGHT',10000,'MINPEAKDISTANCE',100); % Mimic 1
            [~,sample1_M2]=findpeaks(((A(2,:))) ,'MINPEAKHEIGHT',10000,'MINPEAKDISTANCE',100); % Mimic 2
            [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)),'MINPEAKHEIGHT',0.8,'MINPEAKDISTANCE',60000); % find event, object-on
            %[~,sample1_N]=findpeaks((B-min(B))/ range(B),'MINPEAKHEIGHT',0.8,'MINPEAKDISTANCE',60000); % find event, object-on
            
            
            
            
            if size(sample1_N,2)<size(M,1) %check if events in csv file correspond to events detected
                [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)) ,'MINPEAKHEIGHT',0.7,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            elseif size(sample1_N,2)>size(M,1) %check if events in csv file correspond to events detected
                [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)) ,'MINPEAKHEIGHT',0.9,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
            end
            AUX23=0.5;
            while size(sample1_N,2) ~= size(M,1)
                if size(sample1_N,2)<size(M,1) %check if events in csv file correspond to events detected
                    [~,sample1_N]=findpeaks((A(4,:)-min(A(4,:)))/ range(A(4,:)) ,'MINPEAKHEIGHT',AUX23,'MINPEAKDISTANCE',3524144); % find event, object-on in next channel
                elseif size(sample1_N,2)>size(M,1) %check if events in csv file correspond to events detected
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
   
                elseif (isempty(sample1_M2)==1 || size(sample1_M2,2)<500) && isempty(sample1_M1)==0     % mimic present at mimic 1

                   a=1; Freq1post_1=[]; AUX3=[]; AUX4=[];
                   
                   for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
                    end
                    
                    for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
                    end
                    
                    for j=1:size(sample1_N,2)
                        [idx2,idx]=min(abs(sample1_M1-sample1_N(j)));
                        
                        Freq1post_1(a,1:5)=(sample1_M1(idx:idx+4)-sample1_N(j))./samplerate1;
                        AUX3(a,1)=M(j);
                        AUX4(a,1)=i;
                        Freq1post_1(a,6)=Freq1post_1(a,1)-rand/10;
                        [idx2,idx3]=min(abs(sample1_FM1M2-sample1_M1(idx)));
                        Freq1post_1(a,7)=1/(nanmean(diff(sample1_FM1M2(idx3+3:idx3+13)))/samplerate1);
                        [idx12,idx13]=min(diff(sample1_FM1M2(idx3+3:idx3+13)));
                        Freq1post_1(a,8)=(sample1_FM1M2(idx3+idx13)-sample1_M1(idx))/samplerate1;
                        a=a+1;
                        
                    end
                    Freqtotalpost_1=[Freqtotalpost_1;Freq1post_1];
                    Obj_idx_1=[Obj_idx_1;[AUX3 AUX4]];
                    
                    for k=1:size(sample1_N,2)
                        [idx2,idx3]=min(abs(sample1_FM1M2-sample1_N(k))); %remove mimic from fish + mimic
                        for l=1:20
                            try
                                [idx2,idx]=min(abs(sample1_M1-sample1_FM1M2(idx3-l)));
                                if sample1_M1(idx)>=sample1_FM1M2(idx3-l)
                                    M1(b,1,l)=(sample1_M1(idx-1)-sample1_FM1M2(idx3-l))/samplerate1;
                                    M1(b,2,l)=(sample1_M1(idx)-sample1_FM1M2(idx3-l))/samplerate1;
                                    M1(b,3,l)=(sample1_N(k)-sample1_FM1M2(idx3-l))/samplerate1;
                                else
                                    M1(b,1,l)=(sample1_M1(idx)-sample1_FM1M2(idx3-l))/samplerate1;
                                    M1(b,2,l)=(sample1_M1(idx+1)-sample1_FM1M2(idx3-l))/samplerate1;
                                    M1(b,3,l)=(sample1_N(k)-sample1_FM1M2(idx3-l))/samplerate1;
                                end
                             catch
                                M1(b,1:3,l)=nan;
                            end
                            
                            try
                                [idx2,idx]=min(abs(sample1_M1-sample1_FM1M2(idx3+l)));
                                if sample1_M1(idx)>=sample1_FM1M2(idx3+l)
                                    M1(b,4,l)=(sample1_M1(idx-1)-sample1_FM1M2(idx3+l))/samplerate1;
                                    M1(b,5,l)=(sample1_M1(idx)-sample1_FM1M2(idx3+l))/samplerate1;
                                    M1(b,6,l)=(sample1_N(k)-sample1_FM1M2(idx3+l))/samplerate1;
                                else
                                    M1(b,4,l)=(sample1_M1(idx)-sample1_FM1M2(idx3+l))/samplerate1;
                                    M1(b,5,l)=(sample1_M1(idx+1)-sample1_FM1M2(idx3+l))/samplerate1;
                                    M1(b,6,l)=(sample1_N(k)-sample1_FM1M2(idx3+l))/samplerate1;
                                end
                            catch
                                M1(b,4:6,l)=nan;
                            end

                            
                        end
                        b=b+1;
                    end
                    
    
                elseif (isempty(sample1_M1)==1 || size(sample1_M1,2)<500) && isempty(sample1_M2)==0    % mimic present at mimic 2
                    
                    
                    a=1; Freq1post_2=[]; AUX3=[]; AUX4=[];
                    
                    
                    for k=1:size(sample1_M1,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M1(k)-40 & sample1_FM1M2(:)<=sample1_M1(k)+40)=[];
                    end
                    
                    for k=1:size(sample1_M2,2) %remove mimic from fish + mimic
                        sample1_FM1M2(sample1_FM1M2(:)>=sample1_M2(k)-40 & sample1_FM1M2(:)<=sample1_M2(k)+40)=[];
                    end
                    
                    
                    for j=1:size(sample1_N,2)
                        [idx2,idx]=min(abs(sample1_M2-sample1_N(j)));
                        Freq1post_2(a,1:5)=(sample1_M2(idx:idx+4)-sample1_N(j))./samplerate1;
                        AUX3(a,1)=M(j);
                        AUX4(a,1)=i;
                        Freq1post_2(a,6)=Freq1post_2(a,1)-rand/10;
                        [idx2,idx3]=min(abs(sample1_FM1M2-sample1_M2(idx)));
                        Freq1post_2(a,7)=1/(nanmean(diff(sample1_FM1M2(idx3+3:idx3+13)))/samplerate1);
                        [idx12,idx13]=min(diff(sample1_FM1M2(idx3+3:idx3+13)));
                        Freq1post_2(a,8)=(sample1_FM1M2(idx3+idx13)-sample1_M2(idx))/samplerate1;
                        a=a+1;
                    end
                    Freqtotalpost_2=[Freqtotalpost_2;Freq1post_2];
                    Obj_idx_2=[Obj_idx_2;[AUX3 AUX4]];
           
                    for k=1:size(sample1_N,2)
                        [idx2,idx3]=min(abs(sample1_FM1M2-sample1_N(k))); %remove mimic from fish + mimic
                        for l=1:20
                            try
                                [idx2,idx]=min(abs(sample1_M2-sample1_FM1M2(idx3-l)));
                                if sample1_M2(idx)>=sample1_FM1M2(idx3-l)
                                    M2(c,1,l)=(sample1_M2(idx-1)-sample1_FM1M2(idx3-l))/samplerate1;
                                    M2(c,2,l)=(sample1_M2(idx)-sample1_FM1M2(idx3-l))/samplerate1;
                                    M2(c,3,l)=(sample1_N(k)-sample1_FM1M2(idx3-l))/samplerate1;
                                else
                                    M2(c,1,l)=(sample1_M2(idx)-sample1_FM1M2(idx3-l))/samplerate1;
                                    M2(c,2,l)=(sample1_M2(idx+1)-sample1_FM1M2(idx3-l))/samplerate1;
                                    M2(c,3,l)=(sample1_N(k)-sample1_FM1M2(idx3-l))/samplerate1;
                                end
                            catch
                                M2(c,1:3,l)=nan;
                            end
                            
                            
                            try
                                [idx2,idx]=min(abs(sample1_M2-sample1_FM1M2(idx3+l)));
                                if sample1_M2(idx)>=sample1_FM1M2(idx3+l)
                                    M2(c,4,l)=(sample1_M2(idx-1)-sample1_FM1M2(idx3+l))/samplerate1;
                                    M2(c,5,l)=(sample1_M2(idx)-sample1_FM1M2(idx3+l))/samplerate1;
                                    M2(c,6,l)=(sample1_N(k)-sample1_FM1M2(idx3+l))/samplerate1;
                                else
                                    M2(c,4,l)=(sample1_M2(idx)-sample1_FM1M2(idx3+l))/samplerate1;
                                    M2(c,5,l)=(sample1_M2(idx+1)-sample1_FM1M2(idx3+l))/samplerate1;
                                    M2(c,6,l)=(sample1_N(k)-sample1_FM1M2(idx3+l))/samplerate1;
                                end
                            catch
                                M2(c,4:6,l)=nan;
                            end
                                  
                            
                        end
                        c=c+1; 
                    end   
                   
                end               
            end
        end
    end
end
save([myKsDir(1:end-9),'\Fish_',myKsDir(16),'_',myKsDir(end-7:end),'_MIMICS_data.mat'],'Freqtotalpost_1','Freqtotalpost_2', 'Obj_idx_1','Obj_idx_2','M1','M2')



%%
a=1;  figure;
for i=22:1:53
    AUX1=[]; AUX2=[];
    AUX1(:,:)=M1(Obj_idx_1(:,1)==i,1,:); AUX11=reshape(AUX1,[1,size(AUX1,1)*size(AUX1,2)]);
    AUX2(:,:)=M1(Obj_idx_1(:,1)==i,2,:); AUX12=reshape(AUX2,[1,size(AUX2,1)*size(AUX2,2)]);
    subplot(4,8,a);
    h1=histogram(AUX11,25); h1.Normalization = 'probability'; h1.EdgeAlpha= 0;
    hold on;
    h2=histogram(AUX12,25); h2.Normalization = 'probability'; h2.EdgeAlpha= 0;
    %h1.BinWidth = 0.25;
    
    a=a+1;
    %figure; h = histogram2(AUX11,AUX12,100);
end

a=1;  figure;
for i=22:1:53
    AUX1=[]; AUX2=[];
    AUX1(:,:)=M2(Obj_idx_2(:,1)==i,1,:); AUX11=reshape(AUX1,[1,size(AUX1,1)*size(AUX1,2)]);
    AUX2(:,:)=M2(Obj_idx_2(:,1)==i,2,:); AUX12=reshape(AUX2,[1,size(AUX2,1)*size(AUX2,2)]);
    subplot(4,8,a);
    h1=histogram(AUX11,25); h1.Normalization = 'probability'; h1.EdgeAlpha= 0;
    hold on;
    h2=histogram(AUX12,25); h2.Normalization = 'probability'; h2.EdgeAlpha= 0;
    %h1.BinWidth = 0.25;
    
    a=a+1;
    %figure; h = histogram2(AUX11,AUX12,100);
end

AUX21=[]; AUX22=[];
for j=22:28 
    AUX1=[]; AUX2=[];  
    AUX1(:,:)=Freqtotalpost_1(Obj_idx_1(:,1)==j,[1 6 8])-Freqtotalpost_1(Obj_idx_1(:,1)==j,[6]);
    AUX21=[AUX21;AUX1];  
    AUX2(:,:)=Freqtotalpost_2(Obj_idx_2(:,1)==j,[1 6 8])-Freqtotalpost_2(Obj_idx_2(:,1)==j,[6]);
    AUX22=[AUX22;AUX2];   
end

AUX21(AUX21>=.8)=nan;
figure; boxplot(AUX21,'PlotStyle','compact')

AUX22(AUX22>=.8)=nan;
figure; boxplot(AUX22,'PlotStyle','compact')

figure; plot(AUX22(:,1),AUX22(:,3),'.k');
figure; 
for t=1:size(AUX22,1)
    line([1 2],[AUX22(t,1) AUX22(t,3)]);
    hold on;
end

AUX311=[];
for i=30:1:44
    AUX1=[]; AUX2=[]; AUX111=[]; AUX31=[]; 
    
    AUX1(:,:)=M2(Obj_idx_2(:,1)==i,1,:); 
    AUX111(:,:)=M2(Obj_idx_2(:,1)==i,3,20); 
    AUX31=sum(AUX1>=-0.015&AUX1<=-0.010,2)./AUX111;
    AUX311=[AUX311;[AUX31 Freqtotalpost_2(Obj_idx_2(:,1)==i,7)]];
    
    AUX1=[]; AUX2=[]; AUX111=[]; AUX31=[]; 
    
    AUX1(:,:)=M1(Obj_idx_1(:,1)==i,1,:); 
    AUX111(:,:)=M1(Obj_idx_1(:,1)==i,3,20); 
    AUX31=sum(AUX1>=-0.015&AUX1<=-0.010,2)./AUX111;
    AUX311=[AUX311;[AUX31 Freqtotalpost_1(Obj_idx_1(:,1)==i,7)]];
    
    
end

AUX311(AUX311(:,2)<=20)=nan;
figure; plot(AUX311(:,1),AUX311(:,2),'or')

%%
clearvars;

myKsDir = uigetdir('Z:\locker\Fede\5Fish_new_exp\');
files2=dir([myKsDir, '\*MIMIC*']);
Freqtotalpost_1=[]; Freqtotalpost_2=[]; M1=[]; M2=[]; Obj_idx_1=[]; Obj_idx_2=[];
for i=1:size(files2,1)
    load([myKsDir,'\',files2(i).name])
    if i>1
        Freqtotalpost_1=[Freqtotalpost_1;Freqtotalpost_1];
        Freqtotalpost_2=[Freqtotalpost_2;Freqtotalpost_2];
        M1=[M1;M1]; M2=[M2;M2];
        Obj_idx_1=[Obj_idx_1;Obj_idx_1]; Obj_idx_2=[Obj_idx_2;Obj_idx_2];
    end
end

%%

AUX21=[]; 
for i=1:size(M2,1)
    AUX1=[];
    try
        AUX1=abs(M2(i,6,M2(1,5,:)>0.010 & M2(i,5,:)<0.015));
        AUX21=[AUX21;AUX1(1)];
    catch
        AUX21=[AUX21;3];
    end
end


AUX22=[]; AUX23=[];
for j=30:38
     AUX2=[];  AUX3=[];
    AUX2(:,:)=Freqtotalpost_2(Obj_idx_2(:,1)==j,[1 6 8])-Freqtotalpost_2(Obj_idx_2(:,1)==j,[6]);
    AUX22=[AUX22;AUX2]; 
    
    AUX3=AUX21(Obj_idx_2(:,1)==j,1);
    AUX23=[AUX23;AUX3];
end

% AUX21(AUX21>=.8)=nan;
% figure; boxplot(AUX21,'PlotStyle','compact')
% 
% AUX22(AUX22>=.8)=nan;
% figure; boxplot(AUX22,'PlotStyle','compact')

%figure; plot(AUX23,AUX22(:,3),'or'); axis equal; xlim([0 1]); ylim([0 1]); 

AUXfit=[AUX23, AUX22(:,3)]; AUXfit(AUXfit(:,1)>2,:)=[]; AUXfit(AUXfit(:,2)>2,:)=[];

AUXfity=AUXfit(:,2); AUXfitx=AUXfit(:,1); 
%%

AUX21=[]; 
for i=1:size(M1,1)
    AUX1=[];
    try
        AUX1=abs(M1(i,6,M1(1,5,:)>0.010 & M1(i,5,:)<0.015));
        AUX21=[AUX21;AUX1(1)];
    catch
        AUX21=[AUX21;3];
    end
end


AUX22=[]; AUX23=[];
for j=22:28
     AUX2=[];  AUX3=[];
    AUX2(:,:)=Freqtotalpost_1(Obj_idx_1(:,1)==j,[1 6 8])-Freqtotalpost_1(Obj_idx_1(:,1)==j,[6]);
    AUX22=[AUX22;AUX2]; 
    
    AUX3=AUX21(Obj_idx_1(:,1)==j,1);
    AUX23=[AUX23;AUX3];
end

% AUX21(AUX21>=.8)=nan;
% figure; boxplot(AUX21,'PlotStyle','compact')
% 
% AUX22(AUX22>=.8)=nan;
% figure; boxplot(AUX22,'PlotStyle','compact')

%figure; plot(AUX23,AUX22(:,3),'or'); axis equal; xlim([0 1]); ylim([0 1]); 

AUXfit=[AUX23, AUX22(:,3)]; AUXfit(AUXfit(:,1)>2,:)=[]; AUXfit(AUXfit(:,2)>2,:)=[];

AUXfity=AUXfit(:,2); AUXfitx=AUXfit(:,1); 

