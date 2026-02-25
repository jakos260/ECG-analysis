% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands

function GEOM=invInit_withADJ(varargin)

anisotropyRatio=2;
subj='PPD1';
type='_ventricle';
layfile	='nim64.mla';								
bsmfile=[];
group=[];
subgroup =[];
version = [];
baseDir = 'C:\Users\Peter.Damp2\Documents\Data\geometries\';
pp=1;
while pp<=nargin
	if ischar(varargin{pp})
		key=lower(varargin{pp});
		switch key
			case 'subject'
				subj=varargin{pp+1};pp=pp+2;
			case 'group'
				group=varargin{pp+1};pp=pp+2;  
			case 'subgroup'
				subgroup=varargin{pp+1};pp=pp+2;                 
            case 'version'
				version=varargin{pp+1};pp=pp+2;                 
			case 'type'
				type=['_' varargin{pp+1}];pp=pp+2;
			case 'layfile'
				layfile=varargin{pp+1};pp=pp+2;
			case 'anisotropyratio'
				anisotropyRatio=varargin{pp+1};pp=pp+2;
			case 'bsm'
				bsmfile=varargin{pp+1};pp=pp+2;
            case 'basedir'
                baseDir =varargin{pp+1};pp=pp+2;
		end
	end
end
%%

if ~isempty(group)
    if ~isempty(subgroup)
        heartpath=[baseDir group  '\' subj  '\' subgroup '\' subgroup '_' subj version];       
        distpath=[baseDir group  '\' subj  '\' subgroup '\' subj];       
    else
        heartpath=[baseDir group  '\' subj  '\' subj ];    
        distpath = heartpath;
    end  
else
    heartpath=[baseDir subj  '\' subj ];
    distpath = heartpath;
end

if ~exist([heartpath type '.tri'],'file')
	error('subject does not exist')
end


%% heart geomtries atria or ventricle
GEOM=[];
GEOM.heartpath=heartpath;
GEOM.type=type;
GEOM.subject=subj;
[GEOM.VER,GEOM.ITRI]=loadtri([heartpath type '.tri']);
if ~isempty(strfind(type,'atria'))&& exist( [GEOM.heartpath '_atria' '.tri'],'file')
    [GEOM.AVER,GEOM.AITRI]=loadtri([GEOM.heartpath '_atria' '.tri']);
elseif ~isempty(strfind(type,'ventricle')) && exist( [GEOM.heartpath '_ventricle' '.tri'],'file')
    [GEOM.VVER,GEOM.VITRI]=loadtri([GEOM.heartpath '_ventricle' '.tri']);
end

[GEOM.RVER,GEOM.RITRI]=loadtri([heartpath 'rho.tri']);
[GEOM.LVER,GEOM.LITRI]=loadtri([heartpath 'lho.tri']);

if exist([heartpath 'laplacian.mat'],'file')
    load([heartpath 'laplacian.mat']);
    GEOM.LAPL = L *1000;
end



endoVER=[GEOM.RVER; GEOM.LVER];
GEOM.endoVER=zeros(1,length(GEOM.VER));
for i=1:length(GEOM.VER)
	if any(endoVER(:,1)==GEOM.VER(i,1) & endoVER(:,2)==GEOM.VER(i,2) & endoVER(:,3)==GEOM.VER(i,3))
        GEOM.endoVER(i)=1;
	end
end

endoVER=[GEOM.RVER];
GEOM.RendoVER=zeros(1,length(GEOM.VER));
for i=1:length(GEOM.VER)
	if any(endoVER(:,1)==GEOM.VER(i,1) & ...
		   endoVER(:,2)==GEOM.VER(i,2) & ...
		   endoVER(:,3)==GEOM.VER(i,3))
        GEOM.RendoVER(i)=1;
	end
end

GEOM.endoITRI=zeros(length(GEOM.ITRI),1);
for i=1:length(GEOM.ITRI)
    if sum(GEOM.endoVER(GEOM.ITRI(i,:)))==3
        GEOM.endoITRI(i)=1;
    end
end
GEOM.RendoITRI=zeros(length(GEOM.ITRI),1);
for i=1:length(GEOM.ITRI)
    if sum(GEOM.RendoVER(GEOM.ITRI(i,:)))==3
        GEOM.RendoITRI(i)=1;
    end
end



%% ========================= TORSO ========================
[GEOM.TVER,GEOM.TITRI]=loadtri([heartpath,'tor.tri']);
[GEOM.RLVER,GEOM.RLITRI]=loadtri([heartpath 'rlo.tri']);
[GEOM.LLVER,GEOM.LLITRI]=loadtri([heartpath 'llo.tri']);


GEOM.BSM = loadmat( bsmfile );
if isempty(GEOM.BSM)
    GEOM.BSM =  loadmat( bsmfileAlter );
end
if isempty(GEOM.BSM)
    error('no BSM file found');
elseif max(rms(GEOM.BSM)) < 0.1
    GEOM.BSM = GEOM.BSM * 1000;
end
if size(GEOM.BSM,2)== length(GEOM.TVER) 
    GEOM.BSM = GEOM.BSM';
end
if ~isempty(strfind(layfile,'nim')) && size(GEOM.BSM,1) == 64
    GEOM.BSM = zeromean([GEOM.BSM; zeros(size(GEOM.BSM(1,:)))]);    
end

if size(GEOM.BSM,1)~= length(GEOM.TVER) 
    T=intripol(GEOM.TVER,GEOM.TITRI,1:size(GEOM.BSM,1));
    GEOM.BSM = T * GEOM.BSM;
end

GEOM.LAY=loadmat(layfile);


A=loadmat([heartpath type '_edl.mat']);
if ~isempty(strfind(GEOM.type,'atria'))
    GEOM.AMA=zeromean(32*A(1:length(GEOM.TVER),:));
    GEOM.AMAH=zeromean(32*A(size(A,1) - length(GEOM.VER)+1:end,:));
else
    GEOM.AMA=zeromean( 40*A(1:length(GEOM.TVER),:));
    GEOM.AMAH=zeromean(40*A(size(A,1) - length(GEOM.VER)+1:end,:));
end


% select the nodes of the electrodes default is is assuemed that the first
% nodes correspond to the electrode positions



if exist([heartpath 'elnodes.mat'],'file')    
    elNodes = loadmat([heartpath 'elnodes.mat']');
    AMA = GEOM.AMA;
    BSM= GEOM.BSM;
    VER=GEOM.TVER;
    ITRI=GEOM.TITRI;
  
%     changeNodes=sortrows([[1:length(elNodes)]' elNodes],2);
%     changeNodes=([[1:length(elNodes)]' elNodes]);
    changeNodes = elNodes;
    for i = 1 : length(elNodes)
%         i = changeNodes(j,1);
        oldNode = changeNodes(i);
        
        A = AMA(i,:);
        AMA(i,:) = AMA(oldNode,:);
        AMA(oldNode,:) = A;
        
%         B = BSM(i,:);
%         BSM(i,:) = BSM(oldNode,:);          
%         BSM(oldNode,:) = B;
        
        V = VER(i,:);
        VER(i,:) = VER(oldNode,:);
        VER(oldNode,:) = V;
        
        tITRI= ITRI;
        ITRI(ITRI==oldNode) = i;
        ITRI(tITRI==i) = oldNode;
        if any(changeNodes ==i)
            changeNodes(changeNodes==i) = oldNode;
        end       
    end
    GEOM.TVER=VER;
    GEOM.TITRI =ITRI;
    GEOM.AMA = AMA;
    GEOM.BSM = BSM;
end

if strcmp(layfile,'nim64.mla')
	GEOM.wct=[1 2 65];	
elseif strcmp(layfile,'ams65.mla')|| strcmp(layfile,'markams65.mla')
	GEOM.wct=[63 64 65];
else
    GEOM.wct=97:99;
end
% GEOM.elNodes(17:36)=[];
GEOM.AMA = zeromean(GEOM.AMA);
GEOM.BSM = GEOM.BSM;

GEOM.BSM=zeromean(baselinecor(zeromean(GEOM.BSM)));

% T=intripol(GEOM.TVER,GEOM.TITRI,1:size(GEOM.BSM,1));GEOM.BSM=T*GEOM.BSM;
GEOM.distLeads=findDistElecs(GEOM);



if exist([bsmfile(1:end-3) 'act'],'file')
	GEOM.refact=loadmat([bsmfile(1:end-3) 'act']);
end
if exist([bsmfile(1:end-3) 'vm'],'file')
	GEOM.Vm=loadmat([bsmfile(1:end-3) 'vm']);
end


%% ========================= EDL ==========================

if strcmp(layfile,'nim64.mla')
    GEOM.leadname='nijmegen';
    GEOM.v12=[19 26 34 41 48 54 1 2 65];
	GEOM.wct=[1 2 65];	
	GEOM.barrB24=[7 10 12 4 18 24 25 26 27 31 33 34 36 34 47 46  2 45 54 60 58 57 52 61];
	GEOM.luxA=[61  6  5 13 15 17 19 20 22 25 26 27 28 29 30 32 34 35 36 38 41 43 44 53 54 55  2 59 58 66 63 ];% elec 137 ommitted
elseif strcmp(layfile,'ams65.mla')|| strcmp(layfile,'markams65.mla')
    GEOM.leadname='amsterdam';
    GEOM.v12=[12 18 25 31 40 45 63 64 65 ];
	GEOM.luxA=[57 62 61  9  7 10 12 13 15 17 19 20 21 22 23 24 25 26 27 28 32 33 35 44 45 47 42 49 52 53 59 ];	
	GEOM.barrB24=[2  5  7 1 11 16 17 18 19 28 25 26 27 35 39 38 42 43 45 50 52 51 44 57];
	GEOM.wct=[63 64 65];
elseif strcmp(layfile,'arnhem99.mla')
    GEOM.leadname='BSM_(arnhem_99)';    
    GEOM.v12=[20 28 36 45 53 61 97 98 99];
  	GEOM.wct=[97 98 99];
end
if 0 %useWCT
    wct=GEOM.wct;
    GEOM.AMA=GEOM.AMA-ones(size(GEOM.AMA,1),1)*mean(GEOM.AMA(wct,:));
    GEOM.BSM=GEOM.BSM-ones(size(GEOM.BSM,1),1)*mean(GEOM.BSM(wct,:));
end


%% ============= distance matrices geom====================
if exist([distpath type 's.dst3d'],'file')
	GEOM.ADJ=loadmat([distpath type 's.adj3d']);
	GEOM.DIST=loadmat([distpath type 's.dst3d']);
	GEOM.ADJsurf=loadmat([distpath type 's.adj2d']);	
    GEOM.DISTsurf=loadmat([distpath type 's.dst2d']);	
    GEOM.ADJ2W=loadmat([distpath type 's.adjanis']);
	GEOM.DIST2W=loadmat([distpath type 's.dstanis']);
else
	error('no diatnce matrices defined');
end

%%
buur=graphdist(GEOM.ITRI);
D=GEOM.ADJ;D(buur==1)=0;
maxd=max(max(D)); % approximatly the distance between apex and base

if strfind(type,'ventr')
	rendover=GEOM.RendoVER;
	GEOM.Rfreewallver=rendover;
    for i=1:length(GEOM.endoVER)
        if rendover(i)
            GEOM.Rfreewallver(GEOM.DIST(i,:)<20 & GEOM.endoVER==0 )=1;
            if any(GEOM.endoVER==1 & GEOM.RendoVER==0 & GEOM.DIST(i,:)< 35 )||...
                    min(GEOM.DIST(i,GEOM.RendoVER==0 & GEOM.endoVER==1)) < 35
                GEOM.Rfreewallver(i)=0;
            end
        end
    end
    for i=1:length(GEOM.endoVER)
        if GEOM.Rfreewallver(i) && min(GEOM.DIST(i,GEOM.RendoVER==0 & GEOM.endoVER==1)) < 25
            GEOM.Rfreewallver(i) =0;
        end
    end

	GEOM.Rpurkinjever=GEOM.endoVER;
	GEOM.Lpurkinjever=GEOM.endoVER;

	for i=1:length(GEOM.endoVER)
		if any(GEOM.VER(i,1)==GEOM.RVER(:,1) & ...
			   GEOM.VER(i,2)==GEOM.RVER(:,2) & ...
			   GEOM.VER(i,3)==GEOM.RVER(:,3),1)
            GEOM.Lpurkinjever(i)=0;
        else
            GEOM.Rpurkinjever(i)=0;
        end	
			if GEOM.Lpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd)*0.3
				GEOM.Lpurkinjever(i)=0;
			end
			if GEOM.Rpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd)*0.3
				GEOM.Rpurkinjever(i)=0;
			end
% 		end
	end
	if exist(['init' subj],'file')
		eval(['GEOM=init' subj '(GEOM);'])
	end

	GEOM.purkinjever=GEOM.Lpurkinjever;
	GEOM.purkinjever(GEOM.Rpurkinjever==1)=1; 

	%%
	% determine types
	GEOM.part=zeros(size(GEOM.VER,1),1);

	a=find(GEOM.endoVER==0);
	Repi=find(min(GEOM.DIST(GEOM.endoVER==1&GEOM.RendoVER==1,GEOM.endoVER==0))<19);
	Repi=a(Repi);
	a=find(GEOM.endoVER==1&GEOM.RendoVER==1);
	Rsept=find(min(GEOM.DIST(GEOM.endoVER==1&GEOM.RendoVER==0,a))<20);
	Rsept=a(Rsept);
	a=GEOM.RendoVER;a(Rsept)=0;	Rendo=find(a==1);
	a=find(GEOM.endoVER==1&GEOM.RendoVER==0);
	Lsept=find(min(GEOM.DIST(GEOM.endoVER==1&GEOM.RendoVER==1,a))<20);
	Lsept=a(Lsept);
	a=GEOM.endoVER;a(Repi)=1;Lepi=find(a==0);
	a=GEOM.endoVER;a(Rendo)=0;a(Rsept)=0;a(Lsept)=0;Lendo=find(a==1);

	GEOM.part(Lendo)=5;
	GEOM.part(Rendo)=2;
	GEOM.part(Repi)=1;
	GEOM.part(Rsept)=3;
	GEOM.part(Lsept)=4;
	GEOM.part(Lepi)=6;

	GEOM.parts={'Repi';'Rendo';'Rsept';'Lsept';'Lendo';'Lepi'};
	% figure(1);
	% showAtria(GEOM.VER,GEOM.ITRI,GEOM.part);

else
	GEOM.purkinjever=zeros(size(GEOM.VER,1),1);
	GEOM.notchpot=zeros(size(GEOM.VER,1),1);	
	if exist(['init' subj],'file')
		eval(['GEOM=init' subj '(GEOM);'])
	end

end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function leads=findDistElecs(GEOM)

m=mean(GEOM.VER);
T=linetris(GEOM.TVER,GEOM.TITRI,m,m+[1,0,0]);
i1=find(T(:,2)<0);i2=find(T(:,2)>0);
v1=m+T(i1,5)*[1 0 0];
v2=m+T(i2,5)*[1 0 0];
leads=[];

m=v1 + 0.4*(v2-v1);


dV=[GEOM.TVER(:,1)-m(1) GEOM.TVER(:,2)-m(2) GEOM.TVER(:,3)-m(3)];
nv=norm3d(dV);
% showPatch(GEOM.TVER,GEOM.TITRI,nv,'nodes',find(nv==min(nv)))

T=linetris(GEOM.TVER,GEOM.TITRI,m,m+[0,1,0]);
i=GEOM.TITRI(T(1,1),:);
if T(1,3)>.5, v=i(2);elseif T(1,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];
i=GEOM.TITRI(T(2,1),:);
if T(2,3)>.5, v=i(2);elseif T(2,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];

T=linetris(GEOM.TVER,GEOM.TITRI,m,m+[1,0,0]);
i=GEOM.TITRI(T(1,1),:);
if T(1,3)>.5, v=i(2);elseif T(1,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];
i=GEOM.TITRI(T(2,1),:);
if T(2,3)>.5, v=i(2);elseif T(2,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];

T=linetris(GEOM.TVER,GEOM.TITRI,m,m+[1,1,0]);
i=GEOM.TITRI(T(1,1),:);
if T(1,3)>.5, v=i(2);elseif T(1,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];
i=GEOM.TITRI(T(2,1),:);
if T(2,3)>.5, v=i(2);elseif T(2,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];


T=linetris(GEOM.TVER,GEOM.TITRI,m,m+[1,-1,0]);
i=GEOM.TITRI(T(1,1),:);
if T(1,3)>.5, v=i(2);elseif T(1,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];
i=GEOM.TITRI(T(2,1),:);
if T(2,3)>.5, v=i(2);elseif T(2,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];


% ==============================
m=v1 + 0.25*(v2-v1);

T=linetris(GEOM.TVER,GEOM.TITRI,m,m+[0,0,1]);
i=GEOM.TITRI(T(1,1),:);
if T(1,3)>.5, v=i(2);elseif T(1,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];
i=GEOM.TITRI(T(2,1),:);
if T(2,3)>.5, v=i(2);elseif T(2,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];

T=linetris(GEOM.TVER,GEOM.TITRI,m,m+[0,-1,1]);
i=GEOM.TITRI(T(1,1),:);
if T(1,3)>.5, v=i(2);elseif T(1,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];
i=GEOM.TITRI(T(2,1),:);
if T(2,3)>.5, v=i(2);elseif T(2,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];

T=linetris(GEOM.TVER,GEOM.TITRI,m,m+[0,1,1]);
i=GEOM.TITRI(T(1,1),:);
if T(1,3)>.5, v=i(2);elseif T(1,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];
i=GEOM.TITRI(T(2,1),:);
if T(2,3)>.5, v=i(2);elseif T(2,4)>.5,v=i(3);else v=i(1);end
leads=[leads; v];


% showPatch(GEOM.TVER,GEOM.TITRI,nv,'nodes',leads)
