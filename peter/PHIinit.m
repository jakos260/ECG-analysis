function PHIinit(aaaa)


global AS
global ATRIA
global TORSO
global VENTR

% =============================== AS ==================
AS=[];
% colormaps
AS.TALENS=loadmat('tims.mcm');
AS.AVO=loadmat('avopot.mcm');

%========================= ATRIA ========================
ATRIA=[];
TORSO=[];
VENTR=[];

heartpath='c:\ECG_simulation\geometries\ppd1\hearts\40\PPD1_h40_v1_';
torpath=  'c:\ECG_simulation\geometries\ppd1\torso\PPD1_';

depA=loadmat('a01.src');
repA=depA(:,2);depA=depA(:,1);
pSA=loadasci('C:\ECG_simulation\matlab\peter\inverse\invedl_ppd1\cases\domgeom10\PPD1_case10_atria_dis_drep.pS');
depV=loadmat('casev08.src');
repV=depV(:,2);depV=depV(:,1);
pSV=loadasci('C:\ECG_simulation\matlab\peter\inverse\avo case 8\cases\PPD1_case8comp2.pS');

ATRIA.depA=depA;
ATRIA.repA=repA;
ATRIA.pSA=pSA;
VENTR.depV=depV;
VENTR.repV=repV;
VENTR.pSV=pSV;


[ATRIA.VER,ATRIA.ITRI]=loadtri([heartpath 'atria.tri']);
[ATRIA.RVER,ATRIA.RITRI]=loadtri([heartpath 'endo_rcave.tri']);
[ATRIA.LVER,ATRIA.LITRI]=loadtri([heartpath 'endo_lcave.tri']);

endoVER=[ATRIA.RVER; ATRIA.LVER];
ATRIA.endoVER=zeros(1,length(ATRIA.VER));
for i=1:length(ATRIA.VER)
	if ~isempty(find(endoVER(:,1)==ATRIA.VER(i,1) & ...
                     endoVER(:,2)==ATRIA.VER(i,2) & ...
                     endoVER(:,3)==ATRIA.VER(i,3)))
        ATRIA.endoVER(i)=1;
	end
end
ATRIA.endoITRI=zeros(length(ATRIA.ITRI),1);
for i=1:length(ATRIA.ITRI)
    if sum(ATRIA.endoVER(ATRIA.ITRI(i,:)))==3
        ATRIA.endoITRI(i)=1;
    end
end

    
% [ATRIA.RMVER,ATRIA.RMITRI]=loadtri([heartpath 'rim.tri']);
ATRIA.conduction_velocity=0.92; % m/ms

[VENTR.VER,VENTR.ITRI]=loadtri([heartpath 'ventricle.tri']);
VENTR.endoVER=zeros(1,length(VENTR.VER));
for i=1:length(VENTR.VER)
	if ~isempty(find(endoVER(:,1)==VENTR.VER(i,1) & ...
			endoVER(:,2)==VENTR.VER(i,2) & ...
			endoVER(:,3)==VENTR.VER(i,3)))
		VENTR.endoVER(i)=1;
	end
end
VENTR.endoITRI=zeros(length(VENTR.ITRI),1);
for i=1:length(VENTR.ITRI)
    if sum(VENTR.endoVER(VENTR.ITRI(i,:)))==3
        VENTR.endoITRI(i)=1;
    end
end

% VENTR.ADJ=loadmat([heartpath 'ventricle.adj3d']);
% VENTR.DIST=loadmat([heartpath 'ventricle.dis3d']);
VENTR.DIST2D=loadmat([heartpath 'ventricle.dis2d']);
% [temp,VENTR.DIST2D]=graphdist(VENTR.ITRI,VENTR.VER,4);
%========================= TORSO ========================
[TORSO.TVER,TORSO.TITRI]=loadtri([torpath,'tor300.tri']);
[TORSO.RLVER,TORSO.RLITRI]=loadtri([torpath,'RLung_reduced.tri']);
[TORSO.LLVER,TORSO.LLITRI]=loadtri([torpath,'LLung_reduced.tri']);
TORSO.DIS=loadmat([torpath 'tor300_Dist.mat']);

TORSO.BSM=loadmat([torpath 'onebeat01_02.mat']);
TORSO.BSM=zeromean(TORSO.BSM);
% TORSO.wct=[1 2 143];
% TORSO.wct=[141 103 143];
TORSO.wct=[1 2 191];
%========================= EDL ===========================


A=loadmat([heartpath 'atria_edl.mat']);
A=A-ones(size(A,1),1)*mean(A(1:length(TORSO.TVER),:));
V=loadmat([heartpath 'ventricle_edl.mat']);
V=V-ones(size(V,1),1)*mean(V(1:length(TORSO.TVER),:));
A=40*A;
V=40*V;
lA=length(A);
lV=length(V);
disp('complete')

TORSO.AMA_A=A(1:length(TORSO.TVER),:);
TORSO.AMA_V=V(1:length(TORSO.TVER),:);
tl=length(TORSO.TVER);
ATRIA.AMA_LC=A(tl+1:tl+length(ATRIA.LVER),:);
VENTR.AMA_LC=V(tl+1:tl+length(ATRIA.LVER),:);
tl=tl+length(ATRIA.LVER);
ATRIA.AMA_RC=A(tl+1:tl+length(ATRIA.RVER),:);
VENTR.AMA_RC=V(tl+1:tl+length(ATRIA.RVER),:);
tl=tl+length(ATRIA.RVER);
TORSO.AMA_A_RL=A(tl+1:tl+length(TORSO.RLVER),:);
TORSO.AMA_V_RL=V(tl+1:tl+length(TORSO.RLVER),:);	
tl=tl+length(TORSO.RLVER);
TORSO.AMA_A_LL=A(tl+1:tl+length(TORSO.LLVER),:);
TORSO.AMA_V_LL=V(tl+1:tl+length(TORSO.LLVER),:);
tl=tl+length(TORSO.LLVER);
VENTR.AMA_A=A(tl+1:tl+length(ATRIA.VER),:);
ATRIA.AMA_V=V(tl+1:tl+length(VENTR.VER),:);

ATRIA.AMA_A=A(lA-length(ATRIA.VER)+1:end,:);
VENTR.AMA_V=V(lV-length(VENTR.VER)+1:end,:);

[ATRIA.wd,ATRIA.wdVer]=wallthickness(ATRIA.VER,ATRIA.ITRI);
[VENTR.wd,VENTR.wdVer]=wallthickness(VENTR.VER,VENTR.ITRI);

if exist([heartpath 'atria.distw'])
	ATRIA.DISTW=loadmat([heartpath 'atria.distw']);
	ATRIA.ADJW=loadmat([heartpath 'atria.adjw']);
	ATRIA.DIST2W=loadmat([heartpath 'atria.dist2w']);
	ATRIA.ADJ2W=loadmat([heartpath 'atria.adj2w']);
	ATRIA.ADJ=loadmat([heartpath 'atria.adjsec']);
	ATRIA.DIST=loadmat([heartpath 'atria.distsec']);
else
	disp(datestr(clock))
	tic
	[ATRIA.ADJ,ATRIA.DIST]=graphdist(ATRIA.ITRI,ATRIA.VER,3);
	disp(num2str(toc/60))
	disp('updated the distance matrix, now based on second order neigbors')
	disp(datestr(clock))
	tic
    adj=ATRIA.ADJ;
    wd=ATRIA.wd;
    buur=graphdist(ATRIA.ITRI,ATRIA.VER,1);
	disp(num2str(toc/60))
	buur=graphdist(ATRIA.ITRI,ATRIA.VER,1);
	ADJ=ATRIA.ADJ;
	ADJ(buur==0)=2*ADJ(buur==0);
	ATRIA.ADJ2W=ADJ;
	ATRIA.DIST2W=graphdist(ADJ);
	savemat([heartpath 'atria.adj2w'],ATRIA.ADJ2W);
	savemat([heartpath 'atria.dist2w'],ATRIA.DIST2W);
	savemat([heartpath 'atria.adjsec'],ATRIA.ADJ);
	savemat([heartpath 'atria.distsec'],ATRIA.DIST);
	savemat([heartpath 'atria.adjw'],ATRIA.ADJW);
	savemat([heartpath 'atria.distw'],ATRIA.DISTW);
end
if exist([heartpath 'ventr.distw'])
% 	VENTR.DISTW=loadmat([heartpath 'ventr.distw']);
% 	VENTR.ADJW=loadmat([heartpath 'ventr.adjw']);
	VENTR.DIST2W=loadmat([heartpath 'ventr.dist2w']);
	VENTR.ADJ2W=loadmat([heartpath 'ventr.adj2w']);
	VENTR.ADJ=loadmat([heartpath 'ventr.adjsec']);
	VENTR.DIST=loadmat([heartpath 'ventr.distsec']);
else
	disp(datestr(clock))
	tic
	[VENTR.ADJ,VENTR.DIST]=graphdist(VENTR.ITRI,VENTR.VER,3);
	buur=graphdist(VENTR.ITRI,VENTR.VER,1);
	ADJ=VENTR.ADJ;
	ADJ(buur==0)=2*ADJ(buur==0);
% 	for i=1:length(VENTR.endoVER)
% 		a=VENTR.endoVER~=VENTR.endoVER(i);
% 		for j=1:length(a)
% 			if a(j)
% 				ADJ(i,j)=ADJ(i,j)*2;
% 				ADJ(j,i)=ADJ(i,j);
% 			end
% 		end
% 	end
	VENTR.ADJ2W=ADJ;
	VENTR.DIST2W=graphdist(ADJ);
	savemat([heartpath 'ventr.adj2w'],VENTR.ADJ2W);
	savemat([heartpath 'ventr.dist2w'],VENTR.DIST2W);
	savemat([heartpath 'ventr.adjsec'],VENTR.ADJ);
	savemat([heartpath 'ventr.distsec'],VENTR.DIST);
% 	savemat([heartpath 'ventr.adjw'],VENTR.ADJW);
% 	savemat([heartpath 'ventr.distw'],VENTR.DISTW);
% 	
	
end
% [a,b,c,d]=calcActAdj(depV,VENTR.ADJ,1,VENTR.VER,VENTR.ITRI);
% figure(99);clf;colormap(AS.TALENS)
% 	hs=patch('Faces',VENTR.ITRI,'Vertices',VENTR.VER,...
% 					 'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
% 					 'FaceVertexCData',depV,'FaceColor','interp',...
% 					 'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
% 	axis equal off;view(90,0)
% dd=find(d==1);line(VENTR.VER(dd,1),VENTR.VER(dd,2),VENTR.VER(dd,3),'linestyle','none','marker','.')
% 	
%-----------------------------------------------------------------------
%%
function [wd,wopVer]=wallthickness(VER,ITRI)

wd=zeros(length(VER),1);
wopVer=zeros(length(VER),3);
for i=1:length(wd)
	[ti,l]=find(ITRI==i);
	b=cross(VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),2),:),VER(ITRI(ti(:),1),:)-VER(ITRI(ti(:),3),:));
	c=mean(b); c=c./norm(c); % mm
	v1=VER(i,:);
	v2=v1+c;
	TR=linetris(VER,ITRI,v1,v2);
	D=min(TR(TR(:,5)>0.1,5));
	if ~isempty(D) && D < 30
		wd(i)=D;
		wopVer(i,:)=v1+c*D;
	end
end



function c=calcNorm(ver,itri,tri)

if length(tri)==3
	i1=tri(1);
	i2=tri(2);
	i3=tri(3);
	ti=find((itri(:,1)==i1 & itri(:,2)==i2 & itri(:,3)==i3) |...
			(itri(:,1)==i1 & itri(:,2)==i3 & itri(:,3)==i2) |...
			(itri(:,1)==i2 & itri(:,2)==i3 & itri(:,3)==i1) |...
			(itri(:,1)==i2 & itri(:,2)==i1 & itri(:,3)==i3) |...
			(itri(:,1)==i3 & itri(:,2)==i1 & itri(:,3)==i2) |...
			(itri(:,1)==i3 & itri(:,2)==i2 & itri(:,3)==i1));
	c=cross(ver(itri(ti,1),:)-ver(itri(ti,2),:),ver(itri(ti,1),:)-ver(itri(ti,3),:));
	c=c./norm(c); % mm	
elseif length(tri)==2
	i1=tri(1);
	i2=tri(2);
	ti=find((itri(:,1)==i1 & itri(:,2)==i2) |...
			(itri(:,1)==i1 & itri(:,3)==i2) |...
			(itri(:,1)==i2 & itri(:,3)==i1) |...
			(itri(:,1)==i2 & itri(:,2)==i1) |...
			(itri(:,2)==i1 & itri(:,3)==i2) |...
			(itri(:,2)==i2 & itri(:,3)==i1));
	b=cross(ver(itri(ti(:),1),:)-ver(itri(ti(:),2),:),ver(itri(ti(:),1),:)-ver(itri(ti(:),3),:));
	c=mean(b); c=c./norm(c); % mm
else
	[ti,l]=find(itri==tri);
	b=cross(ver(itri(ti(:),1),:)-ver(itri(ti(:),2),:),ver(itri(ti(:),1),:)-ver(itri(ti(:),3),:));
	c=mean(b); c=c./norm(c); % mm
end



	buur=graphdist(VENTR.ITRI,ATRIA.VER,1);
	ADJ=VENTR.ADJ;
	for i=1:length(VENTR.endoVER)
		a=VENTR.endoVER~=VENTR.endoVER(i);
		for j=i+1:length(a)
			if a(j)
				ADJ(i,j)=ADJ(i,j)*2;
				ADJ(j,i)=ADJ(i,j);
			end
		end
	end
	VENTR.ADJ2W=ADJ;
	VENTR.DIST2W=graphdist(ADJ);
	savemat([heartpath 'ventr.adj2w'],VENTR.ADJ2W);
	savemat([heartpath 'ventr.dist2w'],VENTR.DIST2W);

	buur=graphdist(ATRIA.ITRI,ATRIA.VER,1);
	ADJ=ATRIA.ADJ;
	for i=1:length(ATRIA.endoVER)
		a=ATRIA.endoVER~=ATRIA.endoVER(i);
		for j=i+1:length(a)
			if a(j)
				ADJ(i,j)=ADJ(i,j)*2;
				ADJ(j,i)=ADJ(i,j);
			end
		end
	end
	ATRIA.ADJ2W=ADJ;
	ATRIA.DIST2W=graphdist(ADJ);
	savemat([heartpath 'atria.adj2w'],ATRIA.ADJ2W);
	savemat([heartpath 'atria.dist2w'],ATRIA.DIST2W);