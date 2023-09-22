%% read feod for phth and plot

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\8Fish_new_exp\data2\');
files2=dir([myKsDir, '\*TIME_IDX_control*']);
files3=dir([myKsDir, '\*TIME_IDX_test_1*']);
files4=dir([myKsDir, '\*TIME_IDX_test_2*']);
c=0.65;

%%
    ANG=[]; NRTOT=[]; EOD=[];  ME_control=[]; T_control=[]; time_control=[];
    for i=1:size(files2,1)
        load([myKsDir,'\',files2(i).name])
        
        ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   EOD=[EOD EODtime];
        
    end
    
    
    for i=1:size(EOD,2)
        for t=1:32
            AUX=EOD(:,i,t); EODrate1=(diff(AUX)); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end);
            time_control(i,t)=sum(EOD(:,i,t)>=-1 & EOD(:,i,t)<0);
            if sum(isnan(AUX))>=320
                EODrate(1:1801,i,t)=nan;
            else
                try
                    [xs, ys]=FitVal_EI(AUX,EODr1, [-1.3 2.3],0.99999);
                    EODrate(:,i,t)=ys;
                catch
                    EODrate(:,i,t)=nan;
                end
            end
        end
    end
    
    for j=15:19
    [AUX3,AUX4]=max(EODrate(651:1070,:,j));
     AUX4(isnan(AUX3)==1 | AUX3<20)=[];
    time=0:0.002:1;
    MEAUX=time(AUX4); MEAUX(MEAUX<0 | MEAUX>1)=[];
    ME_control=[ME_control, MEAUX];
    
    T_control=[T_control; time_control(time_control(:,j)>2,j)]; 
    end 
    
    
    ANG=[]; NRTOT=[]; EOD=[];  ME_M1=[];  time_M1=[]; T_M1=[];
    for i=1:size(files3,1)
        load([myKsDir,'\',files3(i).name])
        ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   EOD=[EOD EODtime];
    end
    
    for i=1:size(EOD,2)
        for t=1:32
            AUX=EOD(:,i,t); EODrate1=(diff(AUX)); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end);
            time_M1(i,t)=sum(EOD(:,i,t)>=-1 & EOD(:,i,t)<0);
            if sum(isnan(AUX))>=320
                EODrate(1:1801,i,t)=nan;
            else
                try
                    [xs, ys]=FitVal_EI(AUX,EODr1, [-1.3 2.3],0.99999);
                    EODrate(:,i,t)=ys;
                catch
                    EODrate(:,i,t)=nan;
                end
            end
        end
    end
    
    for j=15:19
        [AUX3,AUX4]=max(EODrate(651:1070,:,j));
        AUX4(isnan(AUX3)==1 | AUX3<20)=[];
        time=0:0.002:1;
        MEAUX=time(AUX4); MEAUX(MEAUX<0 | MEAUX>1)=[];
        ME_M1=[ME_M1, MEAUX];
        T_M1=[T_M1; time_M1(time_M1(:,j)>2,j)]; 
    end
    
    
        ANG=[]; NRTOT=[]; EOD=[];  ME_M2=[]; time_M2=[]; T_M2=[];
    for i=1:size(files4,1)
        load([myKsDir,'\',files4(i).name])
        ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   EOD=[EOD EODtime];
    end
    for i=1:size(EOD,2)
        for t=1:32
            AUX=EOD(:,i,t); EODrate1=(diff(AUX)); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end);
            time_M2(i,t)=sum(EOD(:,i,t)>=-1 & EOD(:,i,t)<0);
            if sum(isnan(AUX))>=320
                EODrate(1:1801,i,t)=nan;
            else
                try
                    [xs, ys]=FitVal_EI(AUX,EODr1, [-1.3 2.3],0.99999);
                    EODrate(:,i,t)=ys;
                catch
                    EODrate(:,i,t)=nan;
                end
            end
        end
    end
    
    for j=15:19
    [AUX3,AUX4]=max(EODrate(651:1070,:,j));
     AUX4(isnan(AUX3)==1 | AUX3<20)=[];
    time=0:0.002:1;
    MEAUX=time(AUX4); MEAUX(MEAUX<0 | MEAUX>1)=[];
    ME_M2=[ME_M2, MEAUX];
    T_M2=[T_M2; time_M2(time_M2(:,j)>2,j)]; 
    end 
    
    
    figure; subplot(2,3,1); violinplot(ME_control)
    subplot(2,3,2); violinplot(ME_M1)
    subplot(2,3,3); violinplot(ME_M2)
    
    subplot(2,3,4); violinplot(T_control(T_control<=20)); ylim([0 50]);
    subplot(2,3,5); violinplot(T_M1(T_M1<=20)); ylim([0 50]);
    subplot(2,3,6); violinplot(T_M2(T_M2<=20)); ylim([0 50]);
    
 %%   
    
    
    ANG=[]; NRTOT=[]; EOD=[];
    for i=1:size(files3,1)
        load([myKsDir,'\',files3(i).name])
        ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   EOD=[EOD EODtime];
    end
    for i=1:180
        for t=1:32
            AUX=EOD(:,i,t); EODrate1=(diff(AUX)); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end);
            if sum(isnan(AUX))>=320
                EODrate(1:181,i,t)=nan;
            else
                [xs, ys]=FitVal_EI(AUX,EODr1, [-1.3 2.3],0.99999);
                EODrate(:,i,t)=ys;
            end
        end
    end
    
    MEd=nanmedian(EODrate(:,:,j),2);
    MEdN=1-(MEd./nanmean(EODrate(1:20)));
    Mad=nanstd(EODrate(:,:,j),[],2)/sqrt(length(EODrate(:,:,j))); %std(data)/sqrt(length(data));
    
    subplot(1,3,2); [hl, hp]=boundedline(-1.3:0.02:2.3, MEd,Mad,'-k');
    
    
    ANG=[]; NRTOT=[]; EOD=[];
    for i=1:size(files4,1)
        load([myKsDir,'\',files4(i).name])
        ANG=[ANG;Ang]; NRTOT=[NRTOT Time];   EOD=[EOD EODtime];
    end
    for i=1:180
        for t=1:32
            AUX=EOD(:,i,t); EODrate1=(diff(AUX)); EODr1=1./EODrate1; EODr1(2:end+1)=EODr1(1:end);
            if sum(isnan(AUX))>=320
                EODrate(1:181,i,t)=nan;
            else
                [xs, ys]=FitVal_EI(AUX,EODr1, [-1.3 2.3],0.99999);
                EODrate(:,i,t)=ys;
            end
        end
    end
    
    MEd=nanmedian(EODrate(:,:,j),2);
    MEdN=1-(MEd./nanmean(EODrate(1:20)));
    Mad=nanstd(EODrate(:,:,j),[],2)/sqrt(length(EODrate(:,:,j))); %std(data)/sqrt(length(data));
    
    subplot(1,3,3); [hl, hp]=boundedline(-1.3:0.02:2.3, MEd,Mad,'-k');
%%

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\8Fish_new_exp\data2\');
files2=dir([myKsDir, '\*TIME_IDX_control*']);
files3=dir([myKsDir, '\*TIME_IDX_test_1*']);
files4=dir([myKsDir, '\*TIME_IDX_test_2*']);
c=0.65;


load([myKsDir,'\',files2(1).name])

addpath('D:\KIT3');
%clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\8Fish_new_exp\');
files5=dir([myKsDir, '\*_MIMICS_data.mat']);
%files3=dir([myKsDir, '\*TIME_IDX_test_1*']);
%files4=dir([myKsDir, '\*TIME_IDX_test_2*']);
c=0.65;
load([myKsDir,'\',files5(1).name]);