
% Inverse script voor data Karlrsruhe

% paden toevoegen:
cd 'C:\Users\Jeanne\Dropbox\Jeanne inverse\inverseArno\AM'
addpath(genpath('C:\Users\Jeanne\Dropbox\Jeanne inverse\inverseArno\BEM\inverse'))
addpath('C:\Users\Jeanne\Dropbox\Jeanne inverse\inverseArno\Karlsruhe\matlab\plot')

%% Load alles
% load ventricles:
%[heart.VER heart.ITRI] = loadtri('amhar.tri');

% construct adjecency matrices:
%[heart.ADJsurf, heart.DISTsurf, heart.PATH] = graphdist(heart.ITRI,heart.VER,4);
%[heart.ADJ, heart.DIST, heart.PATH3D]       = graphdist(heart.ITRI,heart.VER,3);

% load endocardial surface geometries:
%[heart.LVER,heart.LITRI] = loadtri([model_KR 'LV_endo.tri']); % Deze zijn
%gemaakt met proper_blood_cavities.m
%[heart.RVER,heart.RITRI] = loadtri([model_KR 'RV_endo.tri']);

% construct anisotropic matrices:
% heart.anisotropyRatio    = 2;                                                   % Not the same anisotropy as used in KarlsRuhe [= 2.6]! [is anisotropy for conduction velocity not conductivity!!]
% heart.ADJ2W              = calcAnisADJ(heart,heart.anisotropyRatio);            % include anisotropy in adjecency matrix
% heart.DIST2W             = graphdist(heart.ADJ2W);                              % include anisotropy in distance matrix

load('ADJ.mat');
load('LAY.mat');
% load('BSM_ectopic.mat');
% load('BSM_ectopic_2.mat');
load('BSM_sinus.mat');

% load thorax geometries:
[TVER,TITRI] = loadtri('thorax.tri');

% load A-matrix:
load('AMA_thorax_inhom.mat');

% construct GEOM:
GEOM.heartpath          = pwd;
GEOM.type               = '_ventricles';
GEOM.VER                = heart.VER;
GEOM.ITRI               = heart.ITRI;
GEOM.TVER               = TVER;
GEOM.TITRI              = TITRI;
GEOM.LVER               = heart.LVER;
GEOM.LITRI              = heart.LITRI;
GEOM.RVER               = heart.RVER;
GEOM.RITRI              = heart.RITRI;
GEOM.ADJ                = heart.ADJ;
GEOM.DIST               = heart.DIST;
GEOM.ADJsurf            = heart.ADJsurf;
GEOM.DISTsurf           = heart.DISTsurf;
GEOM.neigh              = graphdist(GEOM.ITRI);                             %
GEOM.LAY                = LAY;                                              % define LAY file
GEOM.nrfoci             = 1;

GEOM.BSM                = BSM.isotrop_mod;                                  % 1 van de 3 ecg's uit 'BSM_ectopic.mat'
GEOM.AMA                = AMA;
GEOM.anisotropyRatio    = 2;
GEOM.ADJ2W              = heart.ADJ;                                                    % include anisotropy in adjecency matrix
GEOM.DIST2W             = heart.DIST;                                                   % include anisotropy in distance matrix
GEOM.DIST25             = GEOM.DIST2W;                                                  % set anisotropy of initial activation same as later activation.
GEOM                    = selectLeads(GEOM,1:size(GEOM.BSM,1),1);                       % zeromean on BSM and AMA
GEOM.SPECS              = prepareECG(GEOM.BSM, GEOM.LAY, GEOM.neigh, 'documsum',1);

%%
% initial estimate
deprepmode = 4; % 4 voor deprep
[measinit.foci,measinit.dep,measinit.outp,corsinit,rdsinit,deps] = multifociscan(GEOM, 'clusters', GEOM.nrfoci, 'scanmode', deprepmode);

%%
mudep = 1e-5; %murep = mudep;
murep=1.5e-4;
ect=1;
if ect == 1
    measinit.rep	=	measinit.dep + 300;          	% just add 300 ms.
else
    measinit.rep	=	initRep(GEOM,measinit.dep);   % determine repolarization times
end

mrd = 0.15;
% deprepmode = 1; % 4 voor deprep

GEOM.subject    = 'Jeanne1';
GEOM.beat       = '1';

meas = inversefunc(GEOM, measinit.dep, measinit.rep, 'maxiter', 100, 'mudep', mudep, 'murep', murep, 'minrd', mrd, 'mode', deprepmode);

%%
if deprepmode == 1
    sig_length = 1:GEOM.SPECS.qrsduration;
elseif deprepmode == 4
    sig_length = 1:GEOM.SPECS.endtwave;
end

PSIA = lowpassma(GEOM.AMA*getSmode(ones(length(GEOM.VER),1)*(sig_length),meas.depfinal,meas.repfinal,GEOM.SPECS,deprepmode),1);

% construct BSM plot:
h90 = figure(95);
clf;
sigplot(GEOM.BSM(:,1:size(PSIA,2)),'',GEOM.LAY,1/2.5,'b',1,0);
hold on;
sigplot(PSIA,'',GEOM.LAY,1/2.5,'r',1,0);


%% plotten in qtriplot

qtriplot(GEOM.VER,GEOM.ITRI)
qtriplot(measinit.dep) % initial depolarization estimate
qtriplot(measinit.rep) % initial repolarization estimate
qtriplot(meas.depfinal) % final estimate inverse calculated depolarization
qtriplot(meas.repfinal) % final estimate inverse calculated repolarization
qtriplot('funscale') % scale colors if necessary
% qtriplot('delete')