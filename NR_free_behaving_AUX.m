%% test
% %Irate=MA_test(:,:,2).*(PA_test(:,:,2)./nansum((PA_test(:,:,:)),3))+MA_test(:,:,1).*(PA_test(:,:,1)./nansum((PA_test(:,:,:)),3));
%  
% % AUX23(:,:,1)=MA_test(:,:,2).*(PA_test(:,:,2)./nansum((PA_test(:,:,:)),3));
% %  AUX23(:,:,2)=MA_test(:,:,1).*(PA_test(:,:,1)./nansum((PA_test(:,:,:)),3));
% %  Irate=sum(AUX23,3,'omitnan');
% %   Irate((Irate==0))=NaN;
% %  IrateNUM=nansum(PA_test,3);
% 
%   Irate=MA_test(:,:,1);
%   IrateNUM=PA_test(:,:,1);
% 
% x=1:26;y=1:26;
% colormap2 = brewermap(sum(~isnan(reshape(Irate,1,[]))),'Reds');alpha=0.6;alpha2=0.6;
% s=sort(Irate(~isnan(reshape(Irate,1,[]))));circle=10;
% numPoints=100;
% theta=linspace(0,2*pi,numPoints); %100 evenly spaced points between 0 and 2pi
% % imshow(cat(3,Hintergrund,Hintergrund,Hintergrund)); hold on
% figure;
% for i=1:size(x,2)
%     for j=1:size(y,2)
%         if ~isnan(Irate(i,j))
%             col=colormap2((find(Irate(i,j)<=s,1)),:);
%             fPos = get(gcf, 'Position');
%             xl = xlim(); yl = ylim();
% %             w = circle*(xl(2)-xl(1))/fPos(3);
% %             h = circle*(yl(2)-yl(1))/fPos(4);
%             w = circle*IrateNUM(i,j)/fPos(3);
%             h = circle*IrateNUM(i,j)/fPos(4);
%             
%             mx = h*sin(theta); my = h*cos(theta);
%             patch(x(i)+mx*2, y(j)+my*2,col, 'FaceColor', col, 'EdgeColor', 'none');
%         else
%             col=colormap2((find(Irate(i,j)<=s,1)),:);
%             fPos = get(gcf, 'Position');
%             xl = xlim(); yl = ylim();
%             w = 2*(xl(2)-xl(1))/fPos(3);h = 2*(yl(2)-yl(1))/fPos(4);
%             mx = w*sin(theta); my = h*cos(theta);
%             patch(x(i)+mx, y(j)+my*0.8,[0 0 0], 'FaceColor', [0 0 0],  'EdgeColor', [0 0 0]);
%         end
%     end
% end
% hcb=colorbar; axis equal
% set(get(colorbar,'ylabel'),'String', 'EOD rate')
% 
% 
% %% control
% 
% %Irate=MA(:,:,2).*(PA(:,:,2)./nansum((PA(:,:,:)),3))+MA(:,:,1).*(PA(:,:,1)./nansum((PA(:,:,:)),3));
%   
% % AUX23(:,:,1)=MA(:,:,2).*(PA(:,:,2)./nansum((PA(:,:,:)),3));
% %    AUX23(:,:,2)=MA(:,:,1).*(PA(:,:,1)./nansum((PA(:,:,:)),3));
% %    Irate=sum(AUX23,3,'omitnan');
% %     Irate((Irate==0))=NaN;
% %  IrateNUM=nansum(PA,3);
% 
%   Irate=MA(:,:,1);
%   IrateNUM=PA(:,:,1);
% 
% x=1:26;y=1:26;
% colormap2 = brewermap(sum(~isnan(reshape(Irate,1,[]))),'Blues');alpha=0.6;alpha2=0.6;
% s=sort(Irate(~isnan(reshape(Irate,1,[]))));circle=10;
% numPoints=100;
% theta=linspace(0,2*pi,numPoints); %100 evenly spaced points between 0 and 2pi
% % imshow(cat(3,Hintergrund,Hintergrund,Hintergrund)); hold on
% figure;
% for i=1:size(x,2)
%     for j=1:size(y,2)
%         if ~isnan(Irate(i,j))
%             col=colormap2((find(Irate(i,j)<=s,1)),:);
%             fPos = get(gcf, 'Position');
%             xl = xlim(); yl = ylim();
% %             w = circle*(xl(2)-xl(1))/fPos(3);
% %             h = circle*(yl(2)-yl(1))/fPos(4);
% 
%             w = circle*IrateNUM(i,j)/fPos(3);
%             h = circle*IrateNUM(i,j)/fPos(4);
%             
%             mx = h*sin(theta); my = h*cos(theta);
%             patch(x(i)+mx*2, y(j)+my*2,col, 'FaceColor', col, 'EdgeColor', 'none');
%         else
%             col=colormap2((find(Irate(i,j)<=s,1)),:);
%             fPos = get(gcf, 'Position');
%             xl = xlim(); yl = ylim();
%             w = 2*(xl(2)-xl(1))/fPos(3);h = 2*(yl(2)-yl(1))/fPos(4);
%             mx = w*sin(theta); my = h*cos(theta);
%             patch(x(i)+mx, y(j)+my*0.8,[0 0 0], 'FaceColor', [0 0 0],  'EdgeColor', [0 0 0]);
%         end
%     end
% end
% hcb=colorbar; axis equal
% set(get(colorbar,'ylabel'),'String', 'EOD rate')

%% NEW experiment 3 conditions

addpath('D:\KIT3');
clearvars; %close all;
thre=0.15;
%% mimic 1
NRE=[];
myKsDir = uigetdir('Z:\locker\Fede\6Fish_new_exp_2\');
files2=dir([myKsDir, '\*_EOD_data.mat']);

for k=1:size(files2,1)
    load([myKsDir,'\',files2(k).name])
    
    Ang=nan(30,32); NRe=nan(30,32); NRtotal=nan(30,32);
    a=1;
    for i=22:53
        AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
        [AUX1, ~]=find(Obj_idx_1(:,1)==i);
        for j=1:size(AUX1,1)
            NRtotal(j,a)=-mean(FreqTotal_1(AUX1(j),Freqtotalpost_1(AUX1(j),:)<0 & Freqtotalpost_1(AUX1(j),:)>=-1))+mean(FreqTotal_1(AUX1(j),Freqtotalpost_1(AUX1(j),:)>=0 & Freqtotalpost_1(AUX1(j),:)<=1));
            %AUX5=-mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=0 & Freqtotalpost(AUX1(j),:)<=1));
            AUX5=-mean(FreqTotal_1(AUX1(j),Freqtotalpost_1(AUX1(j),:)<0 & Freqtotalpost_1(AUX1(j),:)>=-1))+mean(FreqTotal_1(AUX1(j),Freqtotalpost_1(AUX1(j),:)>=1 & Freqtotalpost_1(AUX1(j),:)<=2));
            
            if NRtotal(j,a)>thre || AUX5>thre
                NRe(j,a)=1;
            else
                NRe(j,a)=0;
            end
        end
        a=a+1;
    end
    
    NRE=[NRE;NRe];
end

POS=[1	26	1	26	1
2	25	2	25	2
3	23	1	23	1
4	22	2	22	2
5	20	1	20	1
6	19	2	19	2
7	17	1	17	1
8	16	2	16	2
9	14	1	14	1
10	13	2	13	2
11	11	1	11	1
12	10	2	10	2
13	8	1	8	1
14	7	2	7	2
15	5	1	5	1
16	4	2	4	2
17	1	4	1	4
18	2	5	2	5
19	1	7	1	7
20	2	8	2	8
21	1	10	1	10
22	2	11	2	11
23	1	13	1	13
24	2	14	2	14
25	1	16	1	16
26	2	17	2	17
27	1	19	1	19
28	2	20	2	20
29	1	22	1	22
30	2	23	2	23
31	1	25	1	25
32	2	26	2	26
];


MA_test1=nan(26,26,1); 
PA_test1=nan(26,26,1);

for j=1:32 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRE(:,j); AUX(isnan(AUX))=[];
    
    MA_test1(POS(j,2),POS(j,3),1)=sum(AUX)/size(AUX,1);
    PA_test1(POS(j,2),POS(j,3),1)=size(AUX,1); 

end


%% mimic 2
NRE=[];
myKsDir = uigetdir('Z:\locker\Fede\5Fish_new_exp\');
files2=dir([myKsDir, '\*_EOD_data.mat']);

for k=1:size(files2,1)
    load([myKsDir,'\',files2(k).name])
    
    Ang=nan(30,32); NRe=nan(30,32); NRtotal=nan(30,32);
    a=1;
    for i=22:53
        AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
        [AUX1, ~]=find(Obj_idx_2(:,1)==i);
        for j=1:size(AUX1,1)
            NRtotal(j,a)=-mean(FreqTotal_2(AUX1(j),Freqtotalpost_2(AUX1(j),:)<0 & Freqtotalpost_2(AUX1(j),:)>=-1))+mean(FreqTotal_2(AUX1(j),Freqtotalpost_2(AUX1(j),:)>=0 & Freqtotalpost_2(AUX1(j),:)<=1));
            %AUX5=-mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=0 & Freqtotalpost(AUX1(j),:)<=1));
            AUX5=-mean(FreqTotal_2(AUX1(j),Freqtotalpost_2(AUX1(j),:)<0 & Freqtotalpost_2(AUX1(j),:)>=-1))+mean(FreqTotal_2(AUX1(j),Freqtotalpost_2(AUX1(j),:)>=1 & Freqtotalpost_2(AUX1(j),:)<=2));
            
            if NRtotal(j,a)>thre || AUX5>thre
                NRe(j,a)=1;
            else
                NRe(j,a)=0;
            end
        end
        a=a+1;
    end
    NRE=[NRE;NRe];
end

POS=[1	26	1	26	1
2	25	2	25	2
3	23	1	23	1
4	22	2	22	2
5	20	1	20	1
6	19	2	19	2
7	17	1	17	1
8	16	2	16	2
9	14	1	14	1
10	13	2	13	2
11	11	1	11	1
12	10	2	10	2
13	8	1	8	1
14	7	2	7	2
15	5	1	5	1
16	4	2	4	2
17	1	4	1	4
18	2	5	2	5
19	1	7	1	7
20	2	8	2	8
21	1	10	1	10
22	2	11	2	11
23	1	13	1	13
24	2	14	2	14
25	1	16	1	16
26	2	17	2	17
27	1	19	1	19
28	2	20	2	20
29	1	22	1	22
30	2	23	2	23
31	1	25	1	25
32	2	26	2	26
];


MA_test2=nan(26,26,1); 
PA_test2=nan(26,26,1);

for j=1:32 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRE(:,j); AUX(isnan(AUX))=[];
    
    MA_test2(POS(j,2),POS(j,3),1)=sum(AUX)/size(AUX,1);
    PA_test2(POS(j,2),POS(j,3),1)=size(AUX,1); 

end

%% control
NRE=[];
myKsDir = uigetdir('Z:\locker\Fede\5Fish_new_exp\');
files2=dir([myKsDir, '\*_EOD_data.mat']);

for k=1:size(files2,1)
    load([myKsDir,'\',files2(k).name])
    
    
    Ang=nan(30,32); NRe=nan(30,32); NRtotal=nan(30,32);
    a=1;
    for i=22:53
        AUX3=[]; AUX2=[]; AUX4=[]; AUX1=[];
        [AUX1, ~]=find(Obj_idx_control(:,1)==i);
        for j=1:size(AUX1,1)
            NRtotal(j,a)=-mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)<0 & Freqtotalpost__control(AUX1(j),:)>=-1))+mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)>=0 & Freqtotalpost__control(AUX1(j),:)<=1));
            %AUX5=-mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)<0 & Freqtotalpost(AUX1(j),:)>=-1))+mean(FreqTotal(AUX1(j),Freqtotalpost(AUX1(j),:)>=0 & Freqtotalpost(AUX1(j),:)<=1));
            AUX5=-mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)<0 & Freqtotalpost__control(AUX1(j),:)>=-1))+mean(FreqTotal_control(AUX1(j),Freqtotalpost__control(AUX1(j),:)>=1 & Freqtotalpost__control(AUX1(j),:)<=2));
            
            if NRtotal(j,a)>thre || AUX5>thre
                NRe(j,a)=1;
            else
                NRe(j,a)=0;
            end
        end
        a=a+1;
    end
    
    NRE=[NRE;NRe];
end

POS=[1	26	1	26	1
2	25	2	25	2
3	23	1	23	1
4	22	2	22	2
5	20	1	20	1
6	19	2	19	2
7	17	1	17	1
8	16	2	16	2
9	14	1	14	1
10	13	2	13	2
11	11	1	11	1
12	10	2	10	2
13	8	1	8	1
14	7	2	7	2
15	5	1	5	1
16	4	2	4	2
17	1	4	1	4
18	2	5	2	5
19	1	7	1	7
20	2	8	2	8
21	1	10	1	10
22	2	11	2	11
23	1	13	1	13
24	2	14	2	14
25	1	16	1	16
26	2	17	2	17
27	1	19	1	19
28	2	20	2	20
29	1	22	1	22
30	2	23	2	23
31	1	25	1	25
32	2	26	2	26
];


MA_control=nan(26,26,1); 
PA_control=nan(26,26,1);

for j=1:32 %obj number
    AUX4=[];  AUX5=[];
    AUX=NRE(:,j); AUX(isnan(AUX))=[];
    
    MA_control(POS(j,2),POS(j,3),1)=sum(AUX)/size(AUX,1);
    PA_control(POS(j,2),POS(j,3),1)=size(AUX,1); 

end

%% plot for test 1 and 2
Irate=MA_test1b(:,:,1);
  IrateNUM=PA_test1(:,:,1);

x=1:26;y=1:26;
colormap2 = brewermap(sum(~isnan(reshape(Irate,1,[]))),'Greys');alpha=0.6;alpha2=0.6;
s=sort(Irate(~isnan(reshape(Irate,1,[]))));circle=10;
numPoints=100;
theta=linspace(0,2*pi,numPoints); %100 evenly spaced points between 0 and 2pi
% imshow(cat(3,Hintergrund,Hintergrund,Hintergrund)); hold on
figure;
for i=1:size(x,2)
    for j=1:size(y,2)
        if ~isnan(Irate(i,j))
            col=colormap2((find(Irate(i,j)<=s,1)),:);
            fPos = get(gcf, 'Position');
            xl = xlim(); yl = ylim();
%             w = circle*(xl(2)-xl(1))/fPos(3);
%             h = circle*(yl(2)-yl(1))/fPos(4);

            w = circle*IrateNUM(i,j)/fPos(3);
            h = circle*IrateNUM(i,j)/fPos(4);
            
            mx = h*sin(theta); my = h*cos(theta);
            patch(x(i)+mx*2, y(j)+my*2,col, 'FaceColor', col, 'EdgeColor', 'none');
        else
            col=colormap2((find(Irate(i,j)<=s,1)),:);
            fPos = get(gcf, 'Position');
            xl = xlim(); yl = ylim();
            w = 2*(xl(2)-xl(1))/fPos(3);h = 2*(yl(2)-yl(1))/fPos(4);
            mx = w*sin(theta); my = h*cos(theta);
            patch(x(i)+mx, y(j)+my*0.8,[0 0 0], 'FaceColor', [0 0 0],  'EdgeColor', [0 0 0]);
        end
    end
end
hcb=colorbar; axis equal
set(get(colorbar,'ylabel'),'String', 'EOD rate')

%% plots for controls

Irate=MA_testcontrolb(:,:,1);
  IrateNUM=PA_control(:,:,1);

x=1:26;y=1:26;
colormap2 = brewermap(sum(~isnan(reshape(Irate,1,[]))),'Greys');alpha=0.6;alpha2=0.6;
s=sort(Irate(~isnan(reshape(Irate,1,[]))));circle=10;
numPoints=100;
theta=linspace(0,2*pi,numPoints); %100 evenly spaced points between 0 and 2pi
% imshow(cat(3,Hintergrund,Hintergrund,Hintergrund)); hold on
figure;
for i=1:size(x,2)
    for j=1:size(y,2)
        if ~isnan(Irate(i,j))
            col=colormap2((find(Irate(i,j)<=s,1)),:);
            fPos = get(gcf, 'Position');
            xl = xlim(); yl = ylim();
%             w = circle*(xl(2)-xl(1))/fPos(3);
%             h = circle*(yl(2)-yl(1))/fPos(4);

            w = circle*IrateNUM(i,j)/fPos(3);
            h = circle*IrateNUM(i,j)/fPos(4);
            
            mx = h*sin(theta); my = h*cos(theta);
            patch(x(i)+mx*2, y(j)+my*2,col, 'FaceColor', col, 'EdgeColor', 'none');
        else
            col=colormap2((find(Irate(i,j)<=s,1)),:);
            fPos = get(gcf, 'Position');
            xl = xlim(); yl = ylim();
            w = 2*(xl(2)-xl(1))/fPos(3);h = 2*(yl(2)-yl(1))/fPos(4);
            mx = w*sin(theta); my = h*cos(theta);
            patch(x(i)+mx, y(j)+my*0.8,[0 0 0], 'FaceColor', [0 0 0],  'EdgeColor', [0 0 0]);
        end
    end
end
hcb=colorbar; axis equal
set(get(colorbar,'ylabel'),'String', 'EOD rate')

%% interpolate

SDS_1=inpaint_nans(MA_test1(:,:,1));
SDS2_1=inpaint_nans(MA_test2(:,:,1));
SDS3_1=inpaint_nans(MA_control(:,:,1));

w     = 5;   % Size of the sliding window (same number of cols and rows in this case)
% Extrapolate values for current window
[Nr,Nc] = size(SDS_1);
Nextra  = 0.5*(w-1);
Ap      = interp2(1:Nc,1:Nr,SDS_1,-Nextra+1:Nc+Nextra,(-Nextra+1:Nr+Nextra).','makima');    % 2D extrapolation must use 'spline' or 'makima' interpolation
% Smooth data with sliding window
H  = ones(w)./w^2;                      % The 2D averaging filter
SDS  = filter2(H,Ap,'valid'); 

Ap      = interp2(1:Nc,1:Nr,SDS2_1,-Nextra+1:Nc+Nextra,(-Nextra+1:Nr+Nextra).','makima');    % 2D extrapolation must use 'spline' or 'makima' interpolation
SDS2  = filter2(H,Ap,'valid'); 
% The smooth resulting matrix

Ap      = interp2(1:Nc,1:Nr,SDS3_1,-Nextra+1:Nc+Nextra,(-Nextra+1:Nr+Nextra).','makima');    % 2D extrapolation must use 'spline' or 'makima' interpolation
SDS3  = filter2(H,Ap,'valid'); 

% figure; 
% ax1=subplot(3,1,1); contourf(xbins(1:end-1),ybins(1:end-1),sds,100,'edgecolor','none'); axis equal; colormap(ax1,brewermap([],'Blues')); %caxis([0 max(max(SDS))]);
% ax2=subplot(3,1,2); contourf(xbins(1:end-1),ybins(1:end-1),sds2,100,'edgecolor','none'); axis equal; colormap(ax2,brewermap([],'Reds')); %caxis([min() max(max(SDS2))]);
% ax3=subplot(3,1,3); contourf(xbins(1:end-1),ybins(1:end-1),sds3,100,'edgecolor','none'); axis equal; colormap(ax3,brewermap([],'Greys')); %caxis([0 max(max(SDS3))]);

% xbins=1:27;
% ybins=1:27;
% figure(100); 
% ax1=subplot(1,3,1); contourf(xbins(1:end-1),ybins(1:end-1),SDS,100,'edgecolor','none'); axis equal; colormap(ax1,brewermap([],'Reds')); %caxis([0 max(max(SDS))]);
% ax2=subplot(1,3,2); contourf(xbins(1:end-1),ybins(1:end-1),SDS2,100,'edgecolor','none'); axis equal; colormap(ax2,brewermap([],'Reds')); %caxis([min() max(max(SDS2))]);
% ax3=subplot(1,3,3); contourf(xbins(1:end-1),ybins(1:end-1),SDS3,100,'edgecolor','none'); axis equal; colormap(ax3,brewermap([],'Blues')); %caxis([0 max(max(SDS3))]);



MA_test1b=nan(26,26); MA_test2b=nan(26,26); MA_testcontrolb=nan(26,26); 
for t=1:size(POS,1)
MA_test1b(POS(t,2),POS(t,3))=SDS(POS(t,2),POS(t,3));
MA_test2b(POS(t,2),POS(t,3))=SDS2(POS(t,2),POS(t,3));
MA_testcontrolb(POS(t,2),POS(t,3))=SDS3(POS(t,2),POS(t,3));
end
