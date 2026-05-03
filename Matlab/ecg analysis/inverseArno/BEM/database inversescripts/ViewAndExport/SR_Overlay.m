
function main
clear all
% close all
display('starting');
h1=[];
h2=[];
zm=0; % zero mean





% load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig03/BasketEnSok/export/varken2443_3BasketEnSok+002_1310_1320_Sinusritme._20111006T135240.mat');
% load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig04/Biosemi/export/varken7358_1_15_SRClosedThorax_4s_64Chan_20111205T113146.mat');
load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig05/Biosemi/export/101930_sr2_10_15_SR_2_ClosedThorax_20111212T143054.mat');

chanarr=1:60;
mla(:,2)=mod(chanarr-1,6)+2;%+1 for the mod trick, only limb leads on 1
mla(:,1)=floor((chanarr-1)/6)+1;
mla(:,3)=chanarr;
mla(31:60,3)=2+mla(31:60,3); % 31 and 32 used as limb leads
mla=[mla;2,1,31;2,8,32;9,1,63;9,8,64];
mla=sortrows(mla,3);
mlas=mla(ismember(mla(:,3),D.SelectChannels),:);
mlas(:,3)=1:size(mlas,1); % renumber channel numbers
mlas(:,1)=mlas(:,1)-min(mlas(:,1))+1;
mlas(:,2)=mlas(:,2)-min(mlas(:,2))+1;
mlas=[max(mlas(:,1:2)),0; mlas];

selbeat=find(DATA.SELBEATS==1,1,'first');
beatclosed=4*DATA.BSMOUT(:,DATA.BEATS(selbeat):DATA.BEATS(selbeat+1));
if zm
    beatclosed=bsxfun(@minus,beatclosed,mean(beatclosed));
end
clear D
clear DATA


% load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig03/BSM/export/varken2443_10_20_Sinus rhythm BSM only Closed thorax_20111006T144047.mat');
% load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig04/Biosemi/export/varken7358_plus2cath_1_15_StartEpiCath_12s_20111205T113024.mat');
% load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig05/Biosemi/export/101930_epimapping_260_290_SrEpiOpenThorax_20111212T144057.mat');
load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig05/Biosemi/export/101930_mapping_1_11_SR_LVEndoMapping_20111212T150048.mat');
% load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig05/Biosemi/export/101930_epimapping_170_180_SR_STartEpiMappingOpenThorax_20111212T154818.mat');
% close all




% for i=1:(length(DATA.BEATS)-1)
figure(1)
set(0,'CurrentFigure',1)
delete(h1);
delete(h2);
i=find(DATA.SELBEATS==1,1,'first');
beatopen=4*DATA.BSMOUT(:,DATA.BEATS(i):DATA.BEATS(i+1));
n=min(size(beatopen,2),size(beatclosed,2));
if zm
    beatopen=bsxfun(@minus,beatopen,mean(beatopen));
end
sig=35;
global PHI2;
%     PHI2=beatopen;
%     sigplotPO(beatclosed,'',mlas,1,'blue');

hold off
sigplot(beatopen(:,1:n),'',mlas,1,'red');
hold on
sigplot(beatclosed(:,1:n),'',mlas,1,'blue');
%     hold off
%     pause(5);


% end



