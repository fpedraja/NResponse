%% read feod for phth and plot

addpath('D:\KIT3');
clearvars; %close all;
myKsDir = uigetdir('Z:\locker\Fede\6Fish_new_exp\data2\');
files2=dir([myKsDir, '\*TIME_IDX_control*']);
files3=dir([myKsDir, '\*TIME_IDX_test_1*']);
files4=dir([myKsDir, '\*TIME_IDX_test_2*']);
c=0.65;

%%
for j=15:19
    figure;
    ANG=[]; NRTOT=[]; EOD=[];
    for i=1:size(files2,1)
        load([myKsDir,'\',files2(i).name])
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
    
    subplot(1,3,1); [hl, hp]=boundedline(-1.3:0.02:2.3, MEd,Mad,'-k');
    
    
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
end