function GEOM = invInit_aj(varargin)
% function GEOM = invInit_aj(varargin)
%
% OPTIONAL INPUT:
% 'subject'         = subject directory [folder] name                        [default = 'PPD1']
% 'group'           = name of group of folders {e.g. 'pigs', 'humans'}       [default = [] ]
% 'type'            = type of surface mesh {e.g. 'ventricles'}               [default = '_ventricle']
% 'layfile'         = layoutfile presentation results BSM                    [default = 'nim65.mla']
% 'anisotropyratio' = anisotropy value                                       [default = 2]
% 'bsm'             = full name of BSM file [*.selecg]                       [default = [] ]
% 'basedir'         = main directory name for Geometries                     [default = '/Users/arnojanssen/Documents/stw/BEM/DATA/Geometries/']
% 'usewct'          = [1] use Wilson Central Terminal specific for 'layfile' [default = 0]
% 'usemean'         = important when BSM files are longer than 4000 samples  [default = 0]
% 'remove50'        = [1] remove 50Hz with method AvO {killhum.m}            [default = 0]
%
% Peter van Dam [PvD]; 2010 november.
% All rights reserved Peacs, Arnhem  the Netherlands
%
% Version version 6-jan-2015, Peter Oosterhoff [PvO] --> Arno Janssen [AJ].
%
% ADD 6-jan-2015: baseDir AJ & spaces in code [for readability]
% Version 2-mar-2015: adapted for AMS PVC STUDY.

% default input:
anisotropyRatio = 2;
subj            = '';
type            = '_ventricle';
layfile         = '';
bsmfile         = [];
group           = [];
baseDir         = '/Users/arnojanssen/Documents/stw/BEM/DATA/Geometries/'; % Arno J's Mac Pro
useWCT          = 0;    % [1] use Wilson's Central Terminal 
usemean         = 0;    % [1] important when BSM files are longer than 4000 samples [SelectBeats.m is still missing!!]
remove50        = 0;    % [1] remove 50Hz from BSM
pp              = 1;

% measured potential difference EDL-layer: 'cardiac potential distributions: A. Van Oosterom'
pot_edl_atr     = 32;
pot_edl_ven     = 40;

% input killhum.m if remove50 = 1: cos+sin-method AvO
f50             = 50;
fsample         = 1000;
mu_ac           = 0.02;

% check input and replace default:
while pp <= nargin,
    if ischar(varargin{pp})
        key = lower(varargin{pp});
        switch key,
            case 'subject'
                subj = varargin{pp+1};              pp = pp+2;
            case 'group'
                group = varargin{pp+1};             pp = pp+2;
            case 'type'
                type = ['_' varargin{pp+1}];        pp = pp+2;
            case 'layfile'
                layfile = varargin{pp+1};           pp = pp+2;
            case 'anisotropyratio'
                anisotropyRatio = varargin{pp+1};   pp = pp+2;
            case 'bsm'
                bsmfile = varargin{pp+1};           pp = pp+2;
            case 'basedir'
                baseDir = varargin{pp+1};           pp = pp+2;
            case 'usewct'
                useWCT = varargin{pp+1};            pp = pp+2;
            case 'usemean'
                usemean = varargin{pp+1};           pp = pp+2;
            case 'remove50'
                remove50 = varargin{pp+1};          pp = pp+2;
        end
    end
end

% select geometry directory for specific subject:
submodel = [subj , '_model'];

if ~isempty(group),
    heartpath       = fullfile(baseDir, group,  subj,  subj );
    anisdistpath    = fullfile(baseDir, group , subj, [subj, num2str(anisotropyRatio), type]);
else
    heartpath       = fullfile(baseDir, submodel , submodel);
    anisdistpath    = fullfile(baseDir , submodel,  [submodel num2str(anisotropyRatio) type ]);
end

if ~exist([heartpath type '.tri'],'file'), error(['subject does not exist:  ' heartpath type '.tri' ]); end

%% heart geomtries atria or ventricle
GEOM                    = [];
GEOM.heartpath          = heartpath;
GEOM.type               = type;
GEOM.subject            = subj;
[GEOM.VER,GEOM.ITRI]    = loadtri([heartpath type '.tri']);
if  exist([heartpath type '.typ'],'file'), GEOM.typ = loadmat([heartpath type '.typ']); end

if strfind(type,'ventricles'),
    [GEOM.VERA,GEOM.ITRIA] = loadtri([heartpath  '_atria.tri']);
    if  exist([heartpath '_atria.typ'],'file'), GEOM.typa = loadmat([heartpath '_atria.typ']); end
end

% area around vertices.
[GEOM.area,tmp] = triareas(GEOM.VER,GEOM.ITRI); 

if ~isempty(strfind(type,'atria'))&& exist( [GEOM.heartpath '_atria' '.tri'],'file')
    [GEOM.AVER,GEOM.AITRI] = loadtri([GEOM.heartpath '_atria' '.tri']);
elseif ~isempty(strfind(type,'ventricle')) && exist( [GEOM.heartpath '_ventricle' '.tri'],'file') % should this be ventricle(s)??
    [GEOM.VVER,GEOM.VITRI] = loadtri([GEOM.heartpath '_ventricle' '.tri']);
end

if exist( [heartpath,'_rcav.tri'],'file')
    [GEOM.RVER,GEOM.RITRI] = loadtri([heartpath '_rcav.tri']);
    [GEOM.LVER,GEOM.LITRI] = loadtri([heartpath '_lcav.tri']);
else
    [GEOM.RVER,GEOM.RITRI] = loadtri([heartpath 'rho.tri']);
    [GEOM.LVER,GEOM.LITRI] = loadtri([heartpath 'lho.tri']);
end

if exist( [heartpath,'_lad.tri'],'file')
    [GEOM.LADVER,GEOM.LADITRI] = loadtri([heartpath '_lad.tri']);
    [GEOM.RCAVER,GEOM.RCAITRI] = loadtri([heartpath '_rca.tri']);
end

if exist([heartpath 'laplacian.mat'],'file')
    load([heartpath 'laplacian.mat']);
    GEOM.LAPL = L*1000;                     % MAGNIFICATION FACTOR!!!! KEEP THIS IN MIND
end

% determine substructers of ventricles, based on typ. 4 classes: endo, Rendo, Lendo, Bendo
if isfield(GEOM,'typ')
    GEOM.endoVER                            = zeros(1,length(GEOM.VER));
    GEOM.RendoVER                           = zeros(1,length(GEOM.VER));
    GEOM.LendoVER                           = zeros(1,length(GEOM.VER));
    GEOM.BendoVER                           = zeros(1,length(GEOM.VER));
    GEOM.endoVER(GEOM.typ==2 | GEOM.typ==3) = 1;
    GEOM.RendoVER(GEOM.typ==3)              = 1;
    GEOM.LendoVER(GEOM.typ==3)              = 1;                        % shouldn't this be 2?!
    GEOM.BendoVER(GEOM.typ>3)               = 1;
else
    endoVER                                 = [GEOM.RVER; GEOM.LVER];
    GEOM.endoVER                            = zeros(1,length(GEOM.VER));
    for i = 1:length(GEOM.VER)
        if any(endoVER(:,1) == GEOM.VER(i,1) & endoVER(:,2) == GEOM.VER(i,2) & endoVER(:,3) == GEOM.VER(i,3))
            GEOM.endoVER(i) = 1;
        end
    end
    
    endoVER         = [GEOM.RVER];
    GEOM.RendoVER   = zeros(1,length(GEOM.VER));
    for i = 1:length(GEOM.VER)
        if any(endoVER(:,1) == GEOM.VER(i,1) & ...
                endoVER(:,2) == GEOM.VER(i,2) & ...
                endoVER(:,3) == GEOM.VER(i,3))
            GEOM.RendoVER(i) = 1;
        end
    end
end

% select tri's complementary to vertices of GEOM.endo
GEOM.endoITRI = zeros(length(GEOM.ITRI),1);
for i = 1:length(GEOM.ITRI)
    if sum(GEOM.endoVER(GEOM.ITRI(i,:))) == 3, GEOM.endoITRI(i) = 1; end
end

% select tri's complementary to vertices of GEOM.Rendo
GEOM.RendoITRI = zeros(length(GEOM.ITRI),1);
for i = 1:length(GEOM.ITRI)
    if sum(GEOM.RendoVER(GEOM.ITRI(i,:))) == 3, GEOM.RendoITRI(i) = 1; end
end

%% ========================= TORSO ========================
if exist( [heartpath,'_thorax.tri'],'file') == 2
    if ~isempty(strfind(layfile,'ams')) &&   exist( [heartpath,'_thoraxAmsterdam.tri'],'file') == 2
        [GEOM.TVER,GEOM.TITRI] = loadtri([heartpath,'_thoraxAmsterdam.tri']);   
    elseif ~isempty(strfind(layfile,'lds')) && exist( [heartpath,'_thorax12_lead.tri'],'file') == 2
        [GEOM.TVER,GEOM.TITRI] = loadtri([heartpath,'_thorax12_lead.tri']);
    else
        [GEOM.TVER,GEOM.TITRI] = loadtri([heartpath,'_thorax.tri']);
        warning('default thorax loaded'); 
    end
    [GEOM.RLVER,GEOM.RLITRI] = loadtri([heartpath '_rlung.tri']);
    [GEOM.LLVER,GEOM.LLITRI] = loadtri([heartpath '_llung.tri']);
else
    [GEOM.TVER,GEOM.TITRI] = loadtri([heartpath,'tor.tri']);
end
[GEOM.Tarea,tmp]    = triareas(GEOM.TVER,GEOM.TITRI);
[GEOM.Tnormal,tmp]  = trinormals(GEOM.TVER,GEOM.TITRI);

GEOM.BSM = loadmat( bsmfile );

if isempty(GEOM.BSM), error('no BSM file found'); end

% remove 50Hz with cos+sin-method AvO
if remove50, GEOM.BSM = killhum(GEOM.BSM,f50,fsample,mu_ac); disp(['removed ' num2str(f50) 'Hz']); end

GEOM.LAY        = loadmat(layfile);                     % visualization grid. dimensions in first row, 1st and 2nd column location in grid, 3rd reference number
GEOM.ECGextra   = GEOM.BSM(max(GEOM.LAY(:,3))+1:end,:); % extra measured channels not represented in visualization grid.

% check if time-signal of bsmfiles is too long.
if size(GEOM.BSM,2) > 4000, error('BSM too long! select beats first.'); end

% loading A-Matrix:
if ~isempty(strfind(layfile,'ams')) &&   exist( [heartpath,'_thoraxAmsterdam.tri'],'file') == 2,
    A = loadmat([heartpath '_thorax_Amsterdam.vedl']);
elseif ~isempty(strfind(layfile,'lds')) && exist( [heartpath,'_thorax12_lead.tri'],'file') == 2,
    GEOM.AMAv12 = pot_edl_ven*loadmat([heartpath '_thorax_12 lead.vedl']);
    A           = loadmat([heartpath type '_v12edl.mat']);
else
    A = loadmat([heartpath type '_edl.mat']);
end

if ~isempty(strfind(GEOM.type,'atria')),
    GEOM.AMA    = zeromean(pot_edl_atr*A(1:length(GEOM.TVER),:));
    GEOM.AMAORG = GEOM.AMA;
    if size(A,1) > length(GEOM.TVER), GEOM.AMAH = zeromean(32*A(size(A,1)-length(GEOM.VER)+1:end,:)); end
else
    GEOM.AMA    = zeromean(pot_edl_ven*A(1:min(size(A,1),length(GEOM.TVER)),:));
    GEOM.AMAORG = GEOM.AMA;
    if size(A,1) > length(GEOM.TVER), GEOM.AMAH = zeromean(pot_edl_ven*A(size(A,1) - length(GEOM.VER)+1:end,:)); end
end

% select the nodes of the electrodes default is assumed that the first
% nodes correspond to the electrode positions

if exist([heartpath 'elnodes.mat'],'file'),
    elNodes     = loadmat([heartpath 'elnodes.mat']');
    AMA         = GEOM.AMA;
    BSM         = GEOM.BSM;
    VER         = GEOM.TVER;
    ITRI        = GEOM.TITRI;
    changeNodes = elNodes;
    for i = 1 : length(elNodes)
        oldNode         = changeNodes(i);
        A               = AMA(i,:);
        AMA(i,:)        = AMA(oldNode,:);
        AMA(oldNode,:)  = A;
        
        B               = BSM(i,:);
        BSM(i,:)        = BSM(oldNode,:);
        BSM(oldNode,:)  = B;
        
        V               = VER(i,:);
        VER(i,:)        = VER(oldNode,:);
        VER(oldNode,:)  = V;
        
        tITRI               = ITRI;
        ITRI(ITRI==oldNode) = i;
        ITRI(tITRI==i)      = oldNode;
        
        if any(changeNodes == i), changeNodes(changeNodes==i) = oldNode; end
    end
    GEOM.TVER   = VER;
    GEOM.TITRI  = ITRI;
    GEOM.AMA    = AMA;
    GEOM.BSM    = BSM;
end

%% ========================= EDL ==========================
GEOM.distLeads = findDistElecs(GEOM);

if strcmp(layfile,'12lds.mla')
    GEOM.leadname   ='12leads';
    GEOM.v12        = 1:9;
    wct             = 1:3;
    GEOM.LAY        = loadmat('9lds.mla');
elseif strcmp(layfile,'ams65.mla')
    GEOM.leadname                   = 'amsterdam';
    GEOM.v12                        = [63 64 65 12 18 25 31 40 45];
    GEOM.v12back                    = [12 18 25 31 40 45 63 64 65 19];      % 12-lead + elec 19 [on back?]
    GEOM.finlay                     = [12 27 9  46 18 57 24 48 22 62 30 15 20 42 59 65 58 47 7  5  38 35 17 61 8  52 32 26 55 13 34 44 21 16 23];
    GEOM.luxA                       = [57 62 61 9  7  10 12 13 15 17 19 20 21 22 23 24 25 26 27 28 32 33 35 44 45 47 42 49 52 53 59 ];
    GEOM.barrB24                    = [2  5  7  1  11 16 17 18 19 28 25 26 27 35 39 38 42 43 45 50 52 51 44 57];
    GEOM.wct                        = [63 64 65];
    GEOM.golem                      = 1:65; 
    GEOM.golem([12 13 14 18 19 46]) = [];
end

% use Wilson's Central Terminal: 
if useWCT,
    wct         = GEOM.wct;
    GEOM.AMA    = GEOM.AMA-ones(size(GEOM.AMA,1),1)*mean(GEOM.AMA(wct,:));
    GEOM.BSM    = GEOM.BSM-ones(size(GEOM.BSM,1),1)*mean(GEOM.BSM(wct,:));
end

%% ============= distance matrices geom ====================
if exist([heartpath type '.dst3d'],'file')                  % load: distance matrices for ventricles. If not present, matrices will be calculated
    GEOM.ADJ        = loadmat([heartpath type '.adj3d']);   % neighbouring vertices in whole ventricles mesh (endo+epi)
    GEOM.DIST       = loadmat([heartpath type '.dst3d']);   % distance vertices in whole ventricles mesh
    GEOM.DISTsurf   = loadmat([heartpath type '.dst2d']);   % distance vertices over surface ventricles mesh
    GEOM.ADJsurf    = loadmat([heartpath type '.adj2d']);   % neighbouring vertices over surface ventricles mesh 
elseif exist([heartpath type '.distsec'],'file')
    GEOM.ADJ        = loadmat([heartpath type '.adjsec']);
    GEOM.DIST       = loadmat([heartpath type '.distsec']);
    GEOM.DISTsurf   = loadmat([heartpath type '.distsurf']);
    GEOM.ADJsurf    = graphdist(GEOM.ITRI,GEOM.VER,4);
else
    disp(datestr(clock))
    tic
    disp('Busy!!!!!')
    [GEOM.ADJ,GEOM.DIST] = graphdist(GEOM.ITRI,GEOM.VER,3);
    disp(num2str(toc/60))
    disp('updated the distance matrix, now based on second order neigbors')
    [GEOM.ADJsurf,GEOM.DISTsurf] = graphdist(GEOM.ITRI,GEOM.VER,4);
    savemat([heartpath type '.adjsec'],GEOM.ADJ);
    savemat([heartpath type '.distsec'],GEOM.DIST);
    savemat([heartpath type '.distsurf'],GEOM.DISTsurf);
end

%% anisotropic matrix calculation
GEOM.anisotropyRatio    = anisotropyRatio;
GEOM.ADJ25              = loadmat([heartpath type '.adjanis']); % including anisotropy ratio of 2.5 in neighbour measure
GEOM.DIST25             = loadmat([heartpath type '.dstanis']); % including anisotropy ratio of 2.5 in distance measure
if anisotropyRatio == 1,
    GEOM.ADJ2W          = GEOM.ADJ;
    GEOM.DIST2W         = GEOM.DIST;
else
    if anisotropyRatio == 2.5 && exist([heartpath type '.adjanis'],'file')
        GEOM.ADJ2W = loadmat([heartpath type '.adjanis']);
    elseif exist([anisdistpath '.adj2w'],'file')
        GEOM.ADJ2W = loadmat([anisdistpath '.adj2w']);
        if length(GEOM.ADJ2W) ~= length(GEOM.ADJ)
            GEOM.ADJ2W = calcAnisADJ(GEOM,anisotropyRatio);
            savemat([anisdistpath '.adj2w'],GEOM.ADJ2W);
        end
    else
        GEOM.ADJ2W = calcAnisADJ(GEOM,anisotropyRatio);
        savemat([anisdistpath '.adj2w'],GEOM.ADJ2W);
    end
    
    if anisotropyRatio == 2.5 && exist([heartpath type '.dstanis'],'file')
        GEOM.DIST2W = loadmat([heartpath type '.dstanis']);    
    elseif exist([anisdistpath '.dist'],'file')
        GEOM.DIST2W = loadmat([anisdistpath '.dist']);
        if length(GEOM.DIST2W) ~= length(GEOM.ADJ)
            GEOM.DIST2W = graphdist(GEOM.ADJ2W);
            savemat([anisdistpath '.dist'],GEOM.DIST2W);
        end
    else
        GEOM.DIST2W = graphdist(GEOM.ADJ2W);
        savemat([anisdistpath '.dist'],GEOM.DIST2W);
    end
end

%%
buur        = graphdist(GEOM.ITRI);
GEOM.neigh  = buur;
D           = GEOM.ADJ;
D(buur==1)  = 0;

if strfind(type,'ventr'),    % determine the right ventricular free wall (outside)
    rendover            = GEOM.RendoVER;
    GEOM.Rfreewallver   = rendover;
    for i = 1:length(GEOM.endoVER)
        if rendover(i),
            GEOM.Rfreewallver(GEOM.DIST(i,:) < 20 & GEOM.endoVER==0) = 1;
            if any(GEOM.endoVER == 1 & GEOM.RendoVER == 0 & GEOM.DIST(i,:) < 15 )|| min(GEOM.DIST(i,GEOM.RendoVER==0 & GEOM.endoVER==1)) < 15,
                GEOM.Rfreewallver(i) = 0;
            end
        end
    end
    
    wd = wallthickness(GEOM.VER,GEOM.ITRI);
    
    for i = 1:length(GEOM.endoVER)
        if GEOM.Rfreewallver(i) && GEOM.RendoVER(i) == 1, % part of free wall and endo
            if  min(GEOM.DIST(i,GEOM.RendoVER==0 & GEOM.endoVER==1)) < 25, GEOM.Rfreewallver(i) = 0; end
        end
    end
    
    if isfield(GEOM,'typ') && ~isempty(GEOM.typ),
        GEOM.Rpostwall  = zeros(size(GEOM.Rfreewallver));
        rvot            = mean(GEOM.VER(GEOM.typ==6,:));
        rvotD           = norm3d(GEOM.VER - ones(length(GEOM.VER),1) * rvot);
        
        ringVer         = mean(GEOM.VER(GEOM.typ==4 ,:));
        dist            = norm3d([GEOM.VER(:,1) - ringVer(1) GEOM.VER(:,2) - ringVer(2) GEOM.VER(:,3) - ringVer(3)]);
        apexi           = find(dist == max( dist( GEOM.typ == 2) ) );
        apexD           = GEOM.DIST(:,apexi);
        
        GEOM.Rpostwall(GEOM.Rfreewallver' == 1 & rvotD > 70 & apexD > 75) = 1;   
    end
       
    GEOM.Rpurkinjever = GEOM.endoVER;
    GEOM.Lpurkinjever = GEOM.endoVER;
    
    for i = 1:length(GEOM.endoVER),
        if any(GEOM.VER(i,1) == GEOM.RVER(:,1) & GEOM.VER(i,2) == GEOM.RVER(:,2) & GEOM.VER(i,3) == GEOM.RVER(:,3),1),
                GEOM.Lpurkinjever(i) = 0;
        else
            GEOM.Rpurkinjever(i) = 0;
        end

        if GEOM.Lpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<= 30, GEOM.Lpurkinjever(i) = 0; end
        if GEOM.Rpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=45,  GEOM.Rpurkinjever(i) = 0; end
    end

    GEOM.purkinjever                        = GEOM.Lpurkinjever;
    GEOM.purkinjever(GEOM.Rpurkinjever==1)  = 1;
    
    %%
    % determine types
    GEOM.part   = zeros(size(GEOM.VER,1),1);
    a           = find(GEOM.endoVER==0);
    Repi        = find(min(GEOM.DIST(GEOM.endoVER==1&GEOM.RendoVER==1,GEOM.endoVER==0))<19);
    Repi        = a(Repi);
    a           = find(GEOM.endoVER==1&GEOM.RendoVER==1);
    Rsept       = find(min(GEOM.DIST(GEOM.endoVER==1&GEOM.RendoVER==0,a))<20);
    Rsept       = a(Rsept);
    a           = GEOM.RendoVER;a(Rsept)=0;	Rendo=find(a==1);
    a           = find(GEOM.endoVER==1&GEOM.RendoVER==0);
    Lsept       = find(min(GEOM.DIST(GEOM.endoVER==1&GEOM.RendoVER==1,a))<20);
    Lsept       = a(Lsept);
    a           = GEOM.endoVER;a(Repi)=1;Lepi=find(a==0);
    a           = GEOM.endoVER;a(Rendo)=0;a(Rsept)=0;a(Lsept)=0;Lendo=find(a==1);
    
    GEOM.part(Lendo)    = 5;
    GEOM.part(Rendo)    = 2;
    GEOM.part(Repi)     = 1;
    GEOM.part(Rsept)    = 3;
    GEOM.part(Lsept)    = 4;
    GEOM.part(Lepi)     = 6;
    GEOM.parts          = {'Repi';'Rendo';'Rsept';'Lsept';'Lendo';'Lepi'};    
else
    GEOM.purkinjever    = zeros(size(GEOM.VER,1),1);
    GEOM.notchpot       = zeros(size(GEOM.VER,1),1);
    if exist(['init' subj],'file'),  eval(['GEOM=init' subj '(GEOM);']); end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function leads = findDistElecs(GEOM) % adapted by AJ 15-jan-2015!

m       = mean(GEOM.VER);                               % determine center of ventricles mesh
T       = linetris(GEOM.TVER,GEOM.TITRI,m,m+[1,0,0]);   % find intersection mesh with line through center in x-direction
i1      = find(T(:,2)<0);i2=find(T(:,2)>0);
v1      = m+T(i1,5)*[1 0 0];                            % calculate intersection line with triangle v1
v2      = m+T(i2,5)*[1 0 0];                            % calculate intersection line with triangle v2
leads   = [];
m040    = v1 + 0.4*(v2-v1);
m025    = v1 + 0.25*(v2-v1);
dV      = [GEOM.TVER(:,1)-m040(1) GEOM.TVER(:,2)-m040(2) GEOM.TVER(:,3)-m040(3)];
nv      = norm3d(dV);
mmat    = [0,1,0; 1,0,0; 1,1,0; 1,-1,0; 0,0,1; 0,-1,1; 0,1,1];

for mvar = 1:7,
    if mvar < 5, m = m040; else m = m025; end
    T = linetris(GEOM.TVER,GEOM.TITRI,m,m+mmat(mvar,:));
    for kvar = 1:2,
        i = GEOM.TITRI(T(kvar,1),:);
        if T(kvar,3)>.5, v = i(2); elseif T(kvar,4)>.5, v = i(3); else v = i(1); end;
        leads   = [leads; v];
    end
end