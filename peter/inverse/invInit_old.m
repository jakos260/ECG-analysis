% Peter van Dam; 2010 november. 
% All rights reserved Peacs, Arnhem  the Netherlands

function GEOM=invInit(varargin)

anisotropyRatio=2;
subj='PPD1';
type='_ventricle';
layfile	='nim65.mla';								
bsmfile=[];
group=[];
baseDir = 'C:\Users\Peter.Damp2\Documents\Data\geometries\';
useWCT = 0;
usemean =0;
remove50 = 0;
pp=1;
while pp<=nargin
	if ischar(varargin{pp})
		key=lower(varargin{pp});
		switch key
			case 'subject'
				subj=varargin{pp+1};pp=pp+2;
			case 'group'
				group=varargin{pp+1};pp=pp+2;                
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
            case 'usewct'
                useWCT = varargin{pp+1};pp=pp+2;
            case 'usemean'
                usemean = varargin{pp+1};pp=pp+2;
            case 'remove50'
                remove50 = varargin{pp+1};pp=pp+2;
		end
	end
end
%%

if ~isempty(group)
    heartpath=[baseDir group  '\' subj  '\' subj ];
    anisdistpath=[baseDir group '\' subj '\' subj num2str(anisotropyRatio) type ];
else
    heartpath=[baseDir subj  '\' subj ];
    anisdistpath=[baseDir '\' subj  '\'  subj num2str(anisotropyRatio) type ];
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

[GEOM.area,tmp] = triareas(GEOM.VER,GEOM.ITRI);
if ~isempty(strfind(type,'atria'))&& exist( [GEOM.heartpath '_atria' '.tri'],'file')
    [GEOM.AVER,GEOM.AITRI]=loadtri([GEOM.heartpath '_atria' '.tri']);
elseif ~isempty(strfind(type,'ventricle')) && exist( [GEOM.heartpath '_ventricle' '.tri'],'file')
    [GEOM.VVER,GEOM.VITRI]=loadtri([GEOM.heartpath '_ventricle' '.tri']);
end

if exist( [heartpath,'_rcav.tri'],'file')
    [GEOM.RVER,GEOM.RITRI]=loadtri([heartpath '_rcav.tri']);
    [GEOM.LVER,GEOM.LITRI]=loadtri([heartpath '_lcav.tri']);
else
    [GEOM.RVER,GEOM.RITRI]=loadtri([heartpath 'rho.tri']);
    [GEOM.LVER,GEOM.LITRI]=loadtri([heartpath 'lho.tri']);   
end

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
if exist( [heartpath,'_thorax.tri'],'file')==2
    if ~isempty(strfind(layfile,'nim')) &&   exist( [heartpath,'_thorax65_Nijmegen.tri'],'file') ==2
        [GEOM.TVER,GEOM.TITRI]=loadtri([heartpath,'_thorax65_Nijmegen.tri']);
    elseif ~isempty(strfind(layfile,'arnhem')) && exist( [heartpath,'_thorax99_Arnhem.tri'],'file') ==2        
        [GEOM.TVER,GEOM.TITRI]=loadtri([heartpath,'_thorax99_Arnhem.tri']);
    else
        [GEOM.TVER,GEOM.TITRI]=loadtri([heartpath,'_thorax.tri']);
    end
    [GEOM.RLVER,GEOM.RLITRI]=loadtri([heartpath '_rlung.tri']);
    [GEOM.LLVER,GEOM.LLITRI]=loadtri([heartpath '_llung.tri']);
else
    [GEOM.TVER,GEOM.TITRI]=loadtri([heartpath,'tor.tri']);
%     [GEOM.RLVER,GEOM.RLITRI]=loadtri([heartpath 'rlo.tri']);
%     [GEOM.LLVER,GEOM.LLITRI]=loadtri([heartpath 'llo.tri']);
end
[GEOM.Tarea,tmp] = triareas(GEOM.TVER,GEOM.TITRI);
[GEOM.Tnormal,tmp] = trinormals(GEOM.TVER,GEOM.TITRI);


GEOM.BSM = loadmat( bsmfile );


if isempty(GEOM.BSM)
    error('no BSM file found');
% elseif max(rms(GEOM.BSM)) < 0.1
%     GEOM.BSM = GEOM.BSM * 1000;
end
if size(GEOM.BSM,2)== length(GEOM.TVER) 
    GEOM.BSM = GEOM.BSM';
end
if remove50
    disp('removed 50Hz')
    GEOM.BSM =killhum(GEOM.BSM,50,1000,0.02);%f50(GEOM.BSM,1000,0.999); %
end

GEOM.LAY=loadmat(layfile);
if size(GEOM.BSM,2) > 2000
    [BSMs,badsigs]=selectBSMBeats(GEOM.BSM(1:size(GEOM.LAY,1)-1,:),GEOM.LAY,'ver',GEOM.TVER,'itri',GEOM.TITRI,'usemean',usemean);     
    GEOM.BSM=baselinecor(BSMs{1});
elseif size(GEOM.BSM,1)~= length(GEOM.TVER) 
    T=intripol(GEOM.TVER,GEOM.TITRI,1:size(GEOM.BSM,1));
    GEOM.BSM = T * GEOM.BSM;
end


A=loadmat([heartpath type '_edl.mat']);

if ~isempty(strfind(GEOM.type,'atria'))
    GEOM.AMAORG = 32*A;
    GEOM.AMA=zeromean(32*A(1:length(GEOM.TVER),:));
    GEOM.AMAH=zeromean(32*A(size(A,1) - length(GEOM.VER)+1:end,:));
else
    GEOM.AMA=zeromean( 40*A(1:length(GEOM.TVER),:));    
    GEOM.AMAORG = 40*A;
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
    changeNodes = elNodes;
    for i = 1 : length(elNodes)
        oldNode = changeNodes(i);
      
        A = AMA(i,:);
        AMA(i,:) = AMA(oldNode,:);
        AMA(oldNode,:) = A;
       
        B = BSM(i,:);
        BSM(i,:) = BSM(oldNode,:);          
        BSM(oldNode,:) = B;
        
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


%% ========================= EDL ==========================
GEOM.distLeads=findDistElecs(GEOM);
if strcmp(layfile,'9lds.mla')
    GEOM.leadname='12leads_Based_nijmegen';
    GEOM.v12=[19 26 34 41 48 54 1 2 65];
elseif strcmp(layfile,'nim64.mla')|| strcmp(layfile,'nim65.mla')
    GEOM.leadname='nijmegen';
    GEOM.v12=[19 26 34 41 48 54 1 2 65];
%     GEOM.v12=[19 26 34 41 47 54 1 2 65];
	GEOM.wct=[1 2 65];	
	GEOM.barrB24=[7 10 12 4 18 24 25 26 27 31 33 34 36 34 47 46  2 45 54 60 58 57 52 61];
	GEOM.luxA=[61  6  5 13 15 17 19 20 22 25 26 27 28 29 30 32 34 35 36 38 41 43 44 53 54 55  2 59 58 66 63 ];% elec 137 ommitted
	GEOM.v12back=[19 26 34 41 48 54 1 2 27];  
	GEOM.finlay=[19	 35	  9	54 26 61 32  51  29   6 40  21 27  2  63  65 62  55 15 17 46  30 25   5  16 58 41 34  59 20  37 53  26 24  23] ;
	if strfind(GEOM.subject,'a'),	GEOM.v12=[GEOM.v12 65];end
elseif strcmp(layfile,'ams65.mla')|| strcmp(layfile,'markams65.mla')
    GEOM.leadname='amsterdam';
    GEOM.v12=[12 18 25 31 40 45  63 64 65 ];
    GEOM.v12back=[12 18 25 31 40 45 63 64 65 19];
% 	finlayorder=[56 106 183 93 73 18 42 155 121 146 75 104 89 12 128 191 96 124 86 24 44 185 41 100 135 47 91 90 127 72 154 60 105  9 152 57    74    92   107   141];
	GEOM.finlay=[12	 27   9 46 18 57 24  48	 22  62 30  15 20 42  59  65 58  47  7  5 38  35 17  61   8 52 32 26  55 13  34 44  21 16  23];
	GEOM.luxA=[57 62 61  9  7 10 12 13 15 17 19 20 21 22 23 24 25 26 27 28 32 33 35 44 45 47 42 49 52 53 59 ];	
	GEOM.barrB24=[2  5  7 1 11 16 17 18 19 28 25 26 27 35 39 38 42 43 45 50 52 51 44 57];
	GEOM.wct=[63 64 65];
	GEOM.golem=1:65;GEOM.golem([12 13 14 18 19 46])=[]; 	
elseif strcmp(layfile,'arnhem99.mla')
    GEOM.leadname='BSM_(arnhem_99)';    
    GEOM.v12=[20 28 36 45 53 61 97 98 99 ];
  	GEOM.wct=[97 98 99];
end

if useWCT
    wct=GEOM.wct;
    GEOM.AMA=GEOM.AMA-ones(size(GEOM.AMA,1),1)*mean(GEOM.AMA(wct,:));
    GEOM.BSM=GEOM.BSM-ones(size(GEOM.BSM,1),1)*mean(GEOM.BSM(wct,:));
end






%% ============= distance matrices geom====================
if exist([heartpath type '.dst3d'],'file')
	GEOM.ADJ=loadmat([heartpath type '.adj3d']);
	GEOM.DIST=loadmat([heartpath type '.dst3d']);
	GEOM.DISTsurf=loadmat([heartpath type '.dst2d']);	
	GEOM.ADJsurf=loadmat([heartpath type '.adj2d']);  

elseif exist([heartpath type '.distsec'],'file')
	GEOM.ADJ=loadmat([heartpath type '.adjsec']);
	GEOM.DIST=loadmat([heartpath type '.distsec']);
	GEOM.DISTsurf=loadmat([heartpath type '.distsurf']);	
	GEOM.ADJsurf=graphdist(GEOM.ITRI,GEOM.VER,4);  
else
	disp(datestr(clock))
	tic
	disp('Busy!!!!!')
	[GEOM.ADJ,GEOM.DIST]=graphdist(GEOM.ITRI,GEOM.VER,3);
	disp(num2str(toc/60))
	disp('updated the distance matrix, now based on second order neigbors')
	[GEOM.ADJsurf,GEOM.DISTsurf]=graphdist(GEOM.ITRI,GEOM.VER,4);
	savemat([heartpath type '.adjsec'],GEOM.ADJ);
	savemat([heartpath type '.distsec'],GEOM.DIST);
	savemat([heartpath type '.distsurf'],GEOM.DISTsurf);
end
%% anisotropic matrix calculation
GEOM.anisotropyRatio = anisotropyRatio;
if anisotropyRatio==1
	GEOM.ADJ2W=GEOM.ADJ;
	GEOM.DIST2W=GEOM.DIST;
else
	if exist([anisdistpath '.adj2w'],'file')
		GEOM.ADJ2W=loadmat([anisdistpath '.adj2w']);
	else
		GEOM.ADJ2W=calcAnisADJ(GEOM,anisotropyRatio);	
		savemat([anisdistpath '.adj2w'],GEOM.ADJ2W);
	end
	if exist([anisdistpath '.dist'],'file')
		GEOM.DIST2W=loadmat([anisdistpath '.dist']);
	else
		GEOM.DIST2W=graphdist(GEOM.ADJ2W);
		savemat([anisdistpath '.dist'],GEOM.DIST2W);
	end
end
%% TODO make ADJ and ADJ2W as sparse as possible

[firstNeighbors,NEIGH] = graphdist(GEOM.ITRI);
GEOM.ADJ2W(GEOM.ADJ>25)=0;
GEOM.ADJ2W(NEIGH >= 2 & GEOM.ADJ./ GEOM.DISTsurf > 0.3 ) = 0;
GEOM.ADJ(GEOM.ADJ>15)=0;
GEOM.ADJ(NEIGH >= 2 & GEOM.ADJ./ GEOM.DISTsurf > 0.3 ) = 0;


% GEOM.ADJ(GEOM.ADJsurf > 0   & firstNeighbors==0) = 0; % remove secondorder neighbors
% GEOM.ADJ2W(GEOM.ADJsurf > 0 & firstNeighbors==0) = 0; % remove secondorder neighbors
% for i=1:length(GEOM.ADJ)
%     minVal = min(GEOM.ADJ(i,firstNeighbors(i,:)==0 & GEOM.ADJ(i,:) > 0));
%     if ~isempty(minVal)
%         GEOM.ADJ(i, GEOM.ADJ(i,:)> minVal & firstNeighbors(i,:)==0) = 0;
%         GEOM.ADJ(GEOM.ADJ(i,:)> minVal & firstNeighbors(i,:)==0,i) = 0;
%     end
%     minVal = min(GEOM.ADJ2W(i,firstNeighbors(i,:)==0 & GEOM.ADJ2W(i,:) > 0));
%     if ~isempty(minVal)
%         GEOM.ADJ2W(i, GEOM.ADJ2W(i,:)> minVal & firstNeighbors(i,:)==0) = 0;
%         GEOM.ADJ2W(GEOM.ADJ2W(i,:)> minVal & firstNeighbors(i,:)==0,i) = 0;
%     end
% end





%%
buur=graphdist(GEOM.ITRI);
D=GEOM.ADJ;D(buur==1)=0;
maxd=max(max(D)); % approximatly the distance between apex and base

if strfind(type,'ventr')
	rendover=GEOM.RendoVER;
	GEOM.Rfreewallver=rendover;
    for i=1:length(GEOM.endoVER)
        if rendover(i)
            GEOM.Rfreewallver(GEOM.DIST(i,:) < 20 & GEOM.endoVER==0 )=1;
            if any(GEOM.endoVER==1 & GEOM.RendoVER==0 & GEOM.DIST(i,:)< 25 )||...
                    min(GEOM.DIST(i,GEOM.RendoVER==0 & GEOM.endoVER==1)) < 25
                GEOM.Rfreewallver(i)=0;
            end
        end
    end
    for i=1:length(GEOM.endoVER)
        if GEOM.Rfreewallver(i) && min(GEOM.DIST(i,GEOM.RendoVER==0 & GEOM.endoVER==1)) < 20
            GEOM.Rfreewallver(i) = 0;
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
		
%         if strfind(GEOM.subject,'go') 
%             if GEOM.Lpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd)*0.25
%                 GEOM.Lpurkinjever(i)=0;
%             end
%             if GEOM.Rpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd)*0.4
% 				GEOM.Rpurkinjever(i)=0;
% 			end
%         elseif strfind(GEOM.subject,'wb') 
%             if GEOM.Lpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd-25)*0.5
%                 GEOM.Lpurkinjever(i)=0;
% 			end
% 			if GEOM.Rfreewallver(i)&& GEOM.Rpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd)*0.64
% 				GEOM.Rpurkinjever(i)=0;
% 			end
% 			if GEOM.Rfreewallver(i)==0 && GEOM.Rpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd)*0.3
% 				GEOM.Rpurkinjever(i)=0;
% 			end
% 		
% 		else
			if GEOM.Lpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<= 30
				GEOM.Lpurkinjever(i)=0;
			end
			if GEOM.Rpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=45
				GEOM.Rpurkinjever(i)=0;
			end
% 		end
	end
% 	if exist(['init' subj],'file')
% 		eval(['GEOM=init' subj '(GEOM);'])
% 	end

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
