%% 
% Inverse functie op amhar.tri met amtor300.tri
cd('C:\Users\Jeanne\Dropbox\Jeanne inverse\inverseArno\AM')
addpath(genpath('C:\Users\Jeanne\Dropbox\Jeanne inverse\inverseArno\BEM\inverse'))
%%

% load ventricles:
[heart.VER, heart.ITRI] = loadtri('amhar.tri');

% load BSM (ECG_PQRST 64 leads)
BSM=loadmat('ampqrst.mes');

load('LAY_64.mat')

% load adjacency matrix, which was made in  construct_adj_amhar.m
load('adj_amhar.mat')

% load A-matrix (.ama converted)
AMA=loadmat('am.ama');

% load thorax geometries:
[TVER,TITRI] = loadtri('amtor300.tri');

%%
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
GEOM.nrfoci             = 3;

GEOM.BSM                = BSM;
GEOM.AMA                = AMA;
GEOM.anisotropyRatio    = heart.anisotropyRatio;
GEOM.ADJ2W              = heart.ADJ2W;                                                    % include anisotropy in adjecency matrix
GEOM.DIST2W             = heart.DIST2W;                                                   % include anisotropy in distance matrix
GEOM.DIST25             = GEOM.DIST2W;                                                  % set anisotropy of initial activation same as later activation.
GEOM                    = selectLeads(GEOM,1:size(GEOM.BSM,1),1);                       % zeromean on BSM and AMA
GEOM.SPECS              = prepareECG(GEOM.BSM, GEOM.LAY, GEOM.neigh, 'documsum',1);


%% initial estimate
deprepmode = 4; % 1 = only depolarisation, 4 = dep + rep, 2=only repolarization - werkt niet helemaal

[measinit.foci,measinit.dep,measinit.outp,corsinit,rdsinit,deps] = multifociscan(GEOM, 'clusters', GEOM.nrfoci, 'scanmode', deprepmode);

%% inverse

mudep = 1e-5; 
murep = mudep;
% murep=1.5e-4;
ect=0;
if ect == 1
    measinit.rep	=	measinit.dep + 300;          	% just add 300 ms.
else
    measinit.rep	=	initRep(GEOM,measinit.dep);   % determine repolarization times
end

mrd = 0.15; % 0.15 = default value

GEOM.subject    = 'Jeanne2';
GEOM.beat       = '1';

meas = inversefunc(GEOM, measinit.dep, measinit.rep, 'maxiter', 100, 'mudep', mudep, 'murep', murep, 'minrd', mrd, 'mode', deprepmode);

%% display results
if deprepmode == 1
    sig_length = 1:GEOM.SPECS.qrsduration;
elseif deprepmode == 4
    sig_length = 1:GEOM.SPECS.endtwave-GEOM.SPECS.onsetqrs;
elseif deprepmode == 2
    sig_length= GEOM.SPECS.onsetqrs+GEOM.SPECS.qrsduration:GEOM.SPECS.endtwave;
end

PSIA = lowpassma(GEOM.AMA*getSmode(ones(length(GEOM.VER),1)*(sig_length),meas.depfinal,meas.repfinal,GEOM.SPECS,deprepmode),1);


% nog define scaling doen!
% sig_length goed gedefinieerd?
% if deprepmode==2
%     maxamp_bsm=max(max(abs(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.onsetqrs+size(PSIA,2)))))/2;
maxamp_bsm=max(max(abs(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.onsetqrs+size(PSIA,2)))))/2;
maxamp_psia=max(max(abs(PSIA)))/2;

% construct BSM plot:
h90 = figure(91);
clf;
sigplot(GEOM.BSM(:,GEOM.SPECS.onsetqrs:GEOM.SPECS.onsetqrs+size(PSIA,2)),'',GEOM.LAY,1/maxamp_bsm,'b',1,0);
hold on;
sigplot(PSIA,'',GEOM.LAY,1/maxamp_psia,'r',1,0);

%% plotten in qtriplot

qtriplot(GEOM.VER,GEOM.ITRI)
qtriplot(measinit.dep) % initial depolarization estimate
qtriplot(measinit.rep) % initial repolarization estimate
qtriplot(meas.depfinal) % final estimate inverse calculated depolarization
qtriplot(meas.repfinal) % final estimate inverse calculated repolarization
qtriplot('funscale') % scale colors if necessary
% qtriplot('delete')
