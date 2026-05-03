clear all
close all
display('starting');
load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig03/BasketEnSok/export/varken2443_3BasketEnSok+002_1310_1320_Sinusritme._20111006T135240.mat');

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
beatopen=4*DATA.BSMOUT(:,DATA.BEATS(selbeat):DATA.BEATS(selbeat+1));
clear D
clear DATA
load('/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig03/BSM/export/varken2443_10_20_Sinus rhythm BSM only Closed thorax_20111006T144047.mat');
selbeat=find(DATA.SELBEATS==1,1,'first');
beatclosed=4*DATA.BSMOUT(:,DATA.BEATS(selbeat):DATA.BEATS(selbeat+1));
n=min(size(beatopen,2),size(beatclosed,2));
beatopen=beatopen(:,1:n);
beatclosed=beatclosed(:,1:n);
figure(1);
sig=35;
global PHI2;
PHI2=beatopen;
sigplotPO(beatclosed,'',mlas,1,'blue');
