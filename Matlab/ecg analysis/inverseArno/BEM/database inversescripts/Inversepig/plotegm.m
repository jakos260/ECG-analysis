clearvars
global measrd measwrd
measwrd=[];
measrd=[];

close all
% geomdir = fullfile('C:\Users\Damp2\Documents\ECG_simulation\STW\Data\geometries\');
geomdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/';

refLeads=1; % 1: zeromean BSM and A-matrix, -2: 9-leads(VR, VL...) apply WCT to A-matrix, BSM already WCT measured. 2: apply WCT to A-matrix and BSM.
casedir		='';
lpass=5%;   % # samples in lowpassma used to filter fw results and BSM!!!

% which leads should be used in the initial esitimate
leadset='all';
group=[];
sub=9; % can be an array now


layfile='pigs_adam.mla';%'pig64.mla';
type='ventricles';
saveCase =1; % 2 = overschrijven ecgspecs, 1 read ecgspecs when available


% INITIALVELOCITY=[0.2,0.4,0.6,0.8,1.0];
% INITIALVELOCITY=[0,0.4,1.0];
initialvelocity=[0.4]%0.4;

maxvelocity=[2.5];%1.6;
velo = 1.0;

% ANIS = [0.01,1,2.5];
anis=0.01;

% INITANIS=[0.01,1,2.5];
initanis=1.0;%2.5;



if sub==2
    subject='ARVCp02';
    %     bsmfile=fullfile('C:\Users\Peter.Damp2\Documents\ECG_simulation\STW\Data\measurements\Pig03\Export\blind\Blind_20111004T193640_File01_beat7.ecg');beat='01';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/ARVC/DATA/ARVCp02/ECG/export/beats/';
    dirout= ['/Users/peteroosterhoff/Documents/Werk/ARVC/DATA/results/' subject '/'];
    refLeads=-2;
elseif sub == 7
    subject='Pig07';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig07/Biosemi/export/AVG/beats/';
    %     dirout= ['.\results\' subject '\'];
    dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
elseif sub == 8
    subject='Pig08';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig08/Biosemi/export/Pacing/AVG/beats/';
    dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
elseif sub == 9
    %         subject='Pig09r';
    subject='Pig09';
    %     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/';
    %     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats3/';BAH,
    %     NOOIT MEER DOEN!!
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/';
    %         bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/SR/AVG/beats/';
    
    %     bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/SRST/beats/';
    dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    %     bsmdir='C:\Users\damp2\Documents\ECG_simulation\STW\Data\measurements\Pig09\ECG\avobeats\';
    %     dirout= ['.\results\' subject '\avobeats\'];
    map2egm=[1:11];% epi1-10, LVEndo1(2-4 not recorded)
elseif sub == 10
    subject='Pig10';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig10/Biosemi/export/AVG/beats/';
    dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    
elseif strcmp(sub,'b01')
    subject='Bucket01';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket01/Biosemi/export/beats/';
    dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    layfile='bucket.lay';
elseif strcmp(sub,'b01s')
    subject='Bucketsmall01';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket01/Biosemi/export/beats/';
    dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    layfile='bucket.lay';
elseif strcmp(sub,'b03')
    subject='Bucket03';
    bsmdir='/Users/peteroosterhoff/Documents/Werk/STW/Data/Bucket03/export/beats/';
    dirout= ['/Users/peteroosterhoff/Documents/Werk/STW/Data/results/' subject '/'];
    layfile='bucket.lay';
    
else
    return
end


% bsmfiles = dir(fullfile(bsmdir, '*.selecg'));
bsmfile='/Users/peteroosterhoff/Documents/Werk/STW/Data/Pig09/Biosemi/export/AVG/beats/lossebeats/pig09_007_580_609_LVLatEpiThrx1SyncVoff_20130211T143819_beat017.selecg';
load('/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/cluster1/FullAVGScan/20140523_im0_wrd0_iV0.4_iAnis1.0Anis0.01/Pig09_pig09_007_580_609_LVLatEpiThrx1SyncVoff_20130211T143819_beat000_ventricles_im0_wrd0_iV0.4_iAnis1.0Anis0.01_23-May-2014init.mat');%init
load('/Users/peteroosterhoff/Documents/Werk/STW/Data/results/Pig09/cluster1/FullAVGScan/20140523_im0_wrd0_iV0.4_iAnis1.0Anis0.01/Pig09_pig09_007_580_609_LVLatEpiThrx1SyncVoff_20130211T143819_beat000_ventricles_im0_wrd0_iV0.4_iAnis1.0Anis0.01_23-May-2014.mat');%inverse
% f=load('');
surftype=loadmat('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/Pig09/Pig09_ventricles.typ');

beat = bsmfile(1:end-7);
[pathstr, name, ext] = fileparts(bsmfile);
orgname = name(1:strfind(name,'beat')-2);
BSM=load([beat(1:strfind(beat,'beats')-1) orgname '.mat']);

ibeat=find(BSM.DATA.SELBEATS,1,'first');
tbeat=BSM.DATA.BEATS(ibeat);


RESULTS =[];

diroutOrg = dirout;

DATE=datestr(now,'yyyymmdd');
% DATE='20140519' % overrule when continuing earlier analysis

% Load all data from selected subject
GEOM=invInit('subject',subject,'layfile',layfile,'bsm',bsmfile,'type',type,'anisotropyRatio',anis,'group',group,'basedir',geomdir);
GEOM.SPECS = prepareECG(GEOM.BSM,GEOM.LAY,'documsum',1,'filename',bsmfile(1:end-7),'dosave',saveCase,'dovstim',1);

t=0:GEOM.SPECS.qrstduration-1;
T=ones(length(GEOM.VER),1)*t;

% PSIA =lowpassma(GEOM.AMAH*getSmode(T,meas.depfinal,meas.repfinal,GEOM.SPECS,4,GEOM),lpass);

[cathpos,realsurf, cathmap]=GetCatheterPosition2(sub,bsmfile);
stimpos=cathpos;
cathsurf=realsurf;
[pnearest, distp, trinearest, focus, distver]=findnearest(stimpos.exactpos,GEOM.VER,GEOM.ITRI,GEOM.typ,cathsurf,1);
[measinitf.foci,measinitf.dep,measinitf.outp] = ...
    multifociscan_publicationPigsPO(GEOM,'clusters',1,'usecor',1,'issinus',0,'initialvelocity',initialvelocity,'maxvelocity',maxvelocity,'focus',focus,'showplots',0,'blmode',0);



PSIAEGM =lowpassma(GEOM.AMAH*getSmode(T,meas.depfinal,meas.repfinal,GEOM.SPECS,4),lpass);
PSIAEGMFOC=lowpassma(GEOM.AMAH*getSmode(T,measinitf.dep,measinitf.rep,GEOM.SPECS,4),lpass);


% PSIAinitnohakkel =GEOM.AMAH*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,4, GEOM);
% PSIAinit =lowpassma(GEOM.AMAH*getSmode(T,measinit.dep,measinit.rep,GEOM.SPECS,4),lpass);

% [cathpos,realsurf, cathmap]=GetCatheterPosition2(sub,bsmfile);
i=1; % late for loop??
vernearest=[];

egmlay=[];
egmlay(:,3)=1:length(map2egm);%cathmap(i).dep);
egmlay(:,2)=egmlay(:,3);
egmlay(:,1)=1;
egmlay=[1 length(map2egm),0 ;egmlay]; % cathmap(i).dep

for iel=1:length(map2egm) %cathmap(i).dep)
    [pnearest,~,trinearest,vernearest(iel),~,la,mu]=findnearest(cathmap(i).exactpos(iel,:),GEOM.VER,GEOM.ITRI,surftype,cathmap(i).surf{iel},1);
end
figure(21);
sscale=.1
% egm=BSM.DATA.EGM(map2egm,tbeat+GEOM.SPECS.onsetqrs:(tbeat+GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint-1));
egm=BSM.DATA.EGM(map2egm,tbeat:(tbeat+GEOM.SPECS.onsetqrs+GEOM.SPECS.time_Jpoint-1));
egm=baselinecor(egm);
sigplot(egm,'',egmlay,sscale,'b');
hold on
% sigplot(PSIA(vernearest,1:GEOM.SPECS.time_Jpoint),'',egmlay,sscale,'r');
sigplot(PSIAEGM(vernearest,1:GEOM.SPECS.time_Jpoint),'',egmlay,sscale,'r');

sigplot(PSIAEGMFOC(vernearest,1:GEOM.SPECS.time_Jpoint),'',egmlay,sscale,'r');








