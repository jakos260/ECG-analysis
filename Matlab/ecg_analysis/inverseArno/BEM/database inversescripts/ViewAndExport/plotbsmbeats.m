disp(['Starting ', mfilename]);
clearvars
selbeatsflag=0; % only plot selected beats
palet='rbkcygmc';
% close all
% filepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/pig07_sr_0_119_Sinusritme_VeelBrom_20120705T115805.mat';
if ~exist('filepath','var')
    [FileName PathName]=uigetfile('*.mat','Select BSM export file','/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/');
    filepath=fullfile(PathName,FileName)
end

load(filepath);
bao=1;
figure(1) 
clf



if selbeatsflag
    DATA.BEATSsel=DATA.BEATS(find(DATA.SELBEATS)');
    nbeats=length(DATA.BEATSsel);
    for i=1:(nbeats)
        hold off
        for j=i:min(nbeats,i+bao-1)
            sigplot(DATA.BSMfilt(:,DATA.BEATSsel(j):(DATA.BEATSsel(j)+1000)),'',DATA.LAY,1,palet(1+mod(j-i,length(palet))),.3,0);
            hold on
        end
        pause
        clf
    end
    
else
    nbeats=length(DATA.BEATS);
    for i=1:bao:(nbeats-1)
        hold off
        for j=i:min(nbeats-1,i+bao-1)
            sigplot(DATA.BSMfilt(:,DATA.BEATS(j):DATA.BEATS(j+1)),'',DATA.LAY,1,palet(1+mod(j-i,length(palet))),.3,0);
            hold on
        end
        pause
        clf
    end
end
