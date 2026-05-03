function GraphOverallSummary
display(['Starting ' mfilename]);
clearvars
close all

% filepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single20131113/OverallSummary.mad';
% filepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single20131113/OverallSummaryEndo.mad';
% filepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/single20131113/OverallSummaryEpi.mad';
if ~exist('filepath','var') || isempty(filepath)
    [filename, pathname]=uigetfile('*.mad','Select OverallSummary file to process','/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/');
    filepath=fullfile(pathname,filename);
end

typestr=filepath(end-7:end-4);
switch typestr
    case 'Endo'
    case 'yEpi'
        typestr='Epi';
    case 'ryLV'
        typestr='LV';
    case 'ryRV'
        typestr='RV';
    otherwise
        typestr='Mixed';
end
typestr=upper(typestr);

load(filepath,'-mat','Savg');
ianis=cell2mat({Savg.anis})==0.1;
initanis=cell2mat({Savg(ianis).initanis});
initialvelocity=cell2mat({Savg(ianis).initialvelocity});
S.initrd=cell2mat({Savg(ianis).Initrd});
S.initcor=cell2mat({Savg(ianis).Initcor});
S.meanmeanabsInitdact=cell2mat({Savg(ianis).meanmeanabsInitdact});
S.maxmeanabsInitdact=cell2mat({Savg(ianis).maxmeanabsInitdact});
S.meandispInitdact=cell2mat({Savg(ianis).meandispInitdact});
S.maxdispInitdact=cell2mat({Savg(ianis).maxdispInitdact});
S.Initdep0=cell2mat({Savg(ianis).Initdep0});
S.initdepstim=cell2mat({Savg(ianis).Initdepstim});
S.initdepstimabs=cell2mat({Savg(ianis).Initdepstimabs});

fields=fieldnames(S);
for i=length(fields):-1:1
    figure(i);
    subplot(1,2,1);
    plotlines(initanis,initialvelocity,getfield(S,fields{i}));
    xlabel('Initial Anisotropy');
    ylabel(fields{i});
    title(sprintf('%s vs Intial Anisotropy for values of intial velocity (%s)',fields{i},typestr));
    subplot(1,2,2);
    plotlines(initialvelocity,initanis,getfield(S,fields{i}));
    xlabel('Initial Velocity');
    ylabel(fields{i});
    title(sprintf('%s vs intial velocity for values of Intial Anisotropy (%s)',fields{i},typestr));
    
end

display(['Finished ' mfilename]);

function plotlines(X,Y,Z)
YV=unique(Y);
cla
for i=1:length(YV)
    ind=find(Y==YV(i));
    plot(X(ind),Z(ind),'-+');
    hold all
    legendlabels{i}=num2str(YV(i));
end
legend(legendlabels{:});

