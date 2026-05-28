disp(['Starting ' mfilename]);
clearvars
[VERv,ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/ajm059/AJM059_ventricles.tri');
[VERt,ITRIt]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/ajm059/AJM059_thorax.tri');
% [VERv,ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/ajm060/AJM060_ventricles.tri');
% [VERt,ITRIt]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/ajm060/AJM060_thorax.tri');
% [VERv,ITRIv]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/ajm069/AJM069_ventricles.tri');
% [VERt,ITRIt]=loadtri('/Users/peteroosterhoff/Documents/Werk/STW/Data/geometries/ajm069/AJM069_thorax.tri');
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm059/FullInverse_mudep1e-1_JopMin';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm059/UsePrev_mudep1e-4_JopMin';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm059/FullInverse_mudep1e-4_JopKnik';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm059/UsePrev_mudep1e-4_Jopknik';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm059/UsePrev_mudep2e-4_JFixed125ms';
basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm059/UsePrev_mudep1e-1_minrd0.17_JFixed125ms';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm060/UsePrev_mudep1e-4';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm060/UsePrev_mudep2e-4';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm069/PrevInit_mudep1e-4';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm069/PrevInit_mudep2e-4';
% basepath='/Users/peteroosterhoff/Documents/Werk/STW/Data/results/ajm069/PrevInit_mudep1e-4_SpecFromBeat1';

basepathbeats='/Users/peteroosterhoff/Documents/Werk/Brugada/DATA/Measurements/AJM059/BSPM/export/beats';
dbeats=dir(fullfile(basepathbeats,'*.selecg'));
dbeatspar=dir(fullfile(basepathbeats,'*.spe'));

signodes=1:65;

qtriplot('delete *');
qtriplot('horizontal 2');
qtriplot('vertical 1');

d=dir(fullfile(basepath,'*2012.mat')); % add 'final' to filename??


T=intripol(VERt,ITRIt,signodes);


depfinal=[];
for i=1:length(d)
    S=load(fullfile(basepath,d(i).name));
    depfinal(i,:)=S.meas.depfinal-min(S.meas.depfinal);
    
    %     depfinal(i,:)=S.meas.depfinal;
end


J=[];
%loop voor beats, signal  V1+ op dummy geometrie
for i=1:length(dbeats)
    bsm=loadmat(fullfile(basepathbeats,dbeats(i).name));
    par=loadmat(fullfile(basepathbeats,dbeatspar(i).name));
    
%     bsm=bsxfun(@minus,bsm,mean(bsm([63 64 65],:),1)); % WTC, check leads
    bsm=baselinecor(bsm);
    VALS= T*bsm(1:length(signodes),:);
    
    lead(i,:)=bsm(17,1:700);
    J(i,:)=VALS(:,par(3));
    
    
end






% qtriplot(VERt,ITRIt,'thorax');
qtriplot(VERv,ITRIv,'Vdep');
% qtriplot(bsxfun(@minus,depfinal',min(bsxfun(@minus,depfinal',depfinal(1,:)'))));% substract diff
qtriplot(depfinal');
qtriplot('funcolor tim');
qtriplot('funscale 0,150');
qtriplot('step 10');
qtriplot('panel Vdep=2,1');


% % plot lead 17, not working
% qtriplot([0 0 0],[],'dummy');
% qtriplot('marker dummy=no');
% % savemat('lead17',lead);
% qtriplot(lead);
% qtriplot('signal dummy=1,0.5,0.3,0.2,0.2 L17)');



% % Difference as left picture
% qtriplot(VERv,ITRIv,'Vdiffdep');
% d=bsxfun(@minus,depfinal',depfinal(1,:)');
% % qtriplot(bsxfun(@minus,d,min(d))); % limit diff to 0
% qtriplot(d);
% % qtriplot('funcolor2 heat');
% qtriplot('panel Vdiffdep=1,1');
% qtriplot('mouse Vdiffdep=fun');


% J-pot on thorax as left picture
qtriplot(VERt/5,ITRIt,'thorax');
qtriplot(J');
qtriplot('fun thorax=2');
qtriplot('funcolor2 heat');
qtriplot('step2 0.05');
qtriplot('funscale2 -1,1');
qtriplot('panel thorax=1,1');
qtriplot('text 0.1,0.8,J-point, 0.1mV step (%column x 10s)');
qtriplot('column 1');

% qtriplot('movie 2');

% qtriplot

disp(['Finished ' mfilename]);