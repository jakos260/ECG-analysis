disp(['Starting ', mfilename]);
clearvars
close all

selbeatsflag=0; % only plot selected beats
palet='rbkcygmc';
% filepath1='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/pig08_1_0_10_SR_Von_ClosedThorax_Inspiration_20121110T205019.mat';
beatstart=[659 1377 2097 2787 4193]; % start QRS



% filepath1='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/pig08_SR_noresp_0_7_SR_Voff_ClosedThorax_20121110T200330.mat';

% filepath2='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/pig08_3_0_35_SR_ThoraxOpen_NoPEEP_Von_20121127T154951.mat';
filepath1='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/pig09_005_283_312_RVApexPostEndoThrx1SyncVoff_20130207T090348_beat000.ecgspecs';
filepath2='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/pig09_005_155_184_RVApexPostEpiThrx1SyncVoff_20130211T133347_beat000.ecgspecs';

% filepath2='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/pig08_3+001_1_10_SR_OpenThoraxAndPericard)Late_20121206T230533.mat';
% filepath2='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/pig08_3+007_1_10_SR_OpenThoraxPericard_Random_20121206T232016.mat';
% filepath2='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/pig08_1_0_10_SR_Von_ClosedThorax_Expiration_20121110T210232.mat';
% close all
% filepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/pig07_sr_0_119_Sinusritme_VeelBrom_20120705T115805.mat';


% filepath1='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/pig07_pace_1050_1109_LV_Apex_Epi_20120627T151225.mat';
% beatstart=[];
% filepath2='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/pig07_pace_1280_1349_LVApexEndo_20120627T160146.mat';


if ~exist('filepath1','var')
    [FileName PathName]=uigetfile('*.ecgspecs','Select  beat  endo','/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/');
    filepath1=fullfile(PathName,FileName);
end

specendo=loadmat(filepath1);

if ~exist('filepath2','var')
    [FileName PathName]=uigetfile('*.ecgspecs','Select BSM export file 2','/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/');
    filepath1=fullfile(PathName,FileName);
end

specepi=loadmat(filepath2);


loadmat


bao=1;
% figure(1)
% clf
lay=loadmat('pig64.mla');
lay=lay(1:63,:);
sl=min(size(S1.DATA.BSMfilt,2),size(S2.DATA.BSMOUT,2))/2;
% sigplot(lowpassma(zeromean(S1.DATA.BSMOUT(1:64,1:sl)),20),'',lay,1,'blue',.3,0);
% hold on
% sigplot(lowpassma(zeromean(S2.DATA.BSMOUT(1:64,1:sl)),20),'',lay,1,'red',.3,0);

bsm1=lowpassma(zeromean(S1.DATA.BSMOUT(1:62,:)),20);
% bsm2=lowpassma(zeromean(S2.DATA.BSMfilt(1:62,:)),20);





figure(2)
r1=rms(bsm1);
plot(r1);
hold on
shift=200;
plot(S1.DATA.BEATS(1:end-1)+shift,r1(S1.DATA.BEATS(1:end-1)+shift),'rx');
% beat1=bsm1(:,4193:4631); % QRS-T


lastbeat=find(S1.DATA.BEATS<sl,1,'last');

% % use selbeats
% S1.DATA.SELBEATS(lastbeat:end)=0; 
beatstart=S1.DATA.BEATS(find(S1.DATA.SELBEATS))+shift; 
beat2i=find(S2.DATA.SELBEATS);
beat2i=beat2i(1);

bsm2=lowpassma(zeromean(S2.DATA.BSMfilt(1:62,S2.DATA.BEATS(beat2i):S2.DATA.BEATS(beat2i+1))),20);

% Use all beats, until lastbeat
% beatstart=S1.DATA.BEATS(1:min(lastbeat-1,1));


% beat1=bsm1(:,4193:4307);
% beat1bl=baselinecor(beat1);
% beatlength=size(beat1,2);
for j=1:length(beatstart)
    
    % beat1=bsm1(:,4193:4307);
    beat1=bsm1(:,beatstart(j):beatstart(j)+300);
    beat1bl=baselinecor(beat1);
    beatlength=size(beat1,2);
    
    
    for i=1:(size(bsm2,2)-beatlength)
        beat2=bsm2(:,i:i+beatlength-1);
        beat2bl=baselinecor(beat2);
        rd(j,i)=norm(beat2bl-beat1bl,'fro')/norm(beat1bl,'fro');
        cc=corrcoef(beat1bl,beat2bl);
        cor(j,i)=cc(1,2);
    end
    [minrd,imin]=min(rd(j,:))
    figure(4)
    clf
    sigplot(beat1bl,'',lay,1,'blue',.3,0);
%     plot(beat1bl(37,:),'blue');
    
    hold on;
    beat2min=baselinecor(bsm2(:,imin:imin+beatlength-1));
    sigplot(beat2min,'',lay,1,'red',.3,0);
%     plot(beat2min(37,:),'red');
    pause
end
figure(3)
plot(rd')
min(rd')
mean(min(rd'))

figure(4)

plot(cor');
min(cor');








% sigplot(S1.DATA.BSMfilt(:,S1.DATA.BEATS(find(S1.DATA.SELBEATS,1)):(S1.DATA.SELBEATS(find(S1.DATA.SELBEATS,1))+1000)),'',S1.DATA.LAY,1,'blue',.3,0);
% hold on
% sigplot(S2.DATA.BSMfilt(:,S2.DATA.SELBEATS(find(S2.DATA.SELBEATS,1)):(S2.DATA.SELBEATS(find(S2.DATA.SELBEATS,1))+1000)),'',S2.DATA.LAY,1,'red',.3,0);



% if selbeatsflag
%     DATA.BEATSsel=DATA.BEATS(find(DATA.SELBEATS)');
%     nbeats=min(length(S1.DATA.BEATSsel);
%     for i=1:(nbeats)
%         hold off
%         for j=i:min(nbeats,i+bao-1)
%             sigplot(DATA.BSMfilt(:,DATA.BEATSsel(j):(DATA.BEATSsel(j)+1000)),'',DATA.LAY,1,palet(1+mod(j-i,length(palet))),.3,0);
%             hold on
%         end
%         pause
%         clf
%     end
%
% else
%     nbeats=length(DATA.BEATS);
%     for i=1:bao:(nbeats-1)
%         hold off
%         for j=i:min(nbeats-1,i+bao-1)
%             sigplot(DATA.BSMfilt(:,DATA.BEATS(j):DATA.BEATS(j+1)),'',DATA.LAY,1,palet(1+mod(j-i,length(palet))),.3,0);
%             hold on
%         end
%         pause
%         clf
%     end
% end
