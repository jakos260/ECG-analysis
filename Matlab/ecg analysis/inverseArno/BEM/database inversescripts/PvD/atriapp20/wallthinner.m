function wallthinner

% merk op dat ik in mm werk. De getallen zijn dan beter leesbaar voor mij.
% Achteraf deel ik weer door 1000

% soothsurface kan waarschijnlijk veel sneller, het gebruikt intripol,
% waarbij maar 1 punt vrijgelaten wodt per slag.


%tris_ra:       1   : 1321           (n= 1321) 
%tris_epi:      1322: 2620           (n= 1299) 
%                                    (n= 2620) 

%tris_ra_endo:  1   :  677           (n=  677)
%tris_ra_epi:    678: 1321           (n=  644)
%                                    (n= 1321)

%tris_la_endo:   1322:2004           (n=  683)
%tris_la_epi:    2005:2620           (n=  616)
%                                    (n= 1299)


inname='atri_pp20_4a';
outname='atri_pp20_4';
[AVER,AITRI]=loadtri([inname '.tri']);
AVER=AVER*1000; % convert to mm

% create epi and endo vertices indices
epiI=[678: 1321 2005:2620];
epi=unique(AITRI(epiI,:));
endoI=[1   :  677 1322:2004];
endo=unique(AITRI(endoI,:));
i=1;
rims=[];
while i~=length(endo)
	if ~isempty(find(epi==endo(i)))
		rims=[rims  endo(i)];
		endo(i)=[];
	else
		i=i+1;
	end
end
i=1;
while i~=length(epi)
	if ~isempty(find(epi(i)==endo))
		epi(i)=[];
	else
		i=i+1;
	end
end

orgAVER=AVER;
nds=[94 82 83 94 82 83];
for i=1:length(nds)
	[a,b]=find(AITRI==nds(i)); 
	neigh=[unique(AITRI(a,:))']; neigh=neigh(neigh~=nds(i));
	AVER(nds(i),:)=mean(AVER(neigh,:));
end
AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry
orgAVER=AVER;


orgAVER=AVER;
disp('smooth heart');
AVER=smoothSurfaceVer(AVER,AITRI,epi);
AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry
orgAVER=AVER;


orgAVER=AVER;
nds=[94 87];
for i=1:length(nds)
	[a,b]=find(AITRI==nds(i)); 
	neigh=[unique(AITRI(a,:))']; neigh=neigh(neigh~=nds(i));
	AVER(nds(i),:)=mean(AVER(neigh,:));
end
AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry
orgAVER=AVER;
TRIstatistics(AVER,AITRI,1000);
% pause
% hier heb ik geprobeerd vooral het linker artoor gevuld te krijgen door
% die grote aftsand langzaam te overbruggen. Ik geloof dat je begint met
% 3.5 cm.

orgAVER=AVER;
for i=30:-1:0
	AVER=thinnerwall(AVER,AITRI,endo,0.935,i+15,2.5,i+11);
	AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry
	orgAVER=AVER;	
end
% Hier wordt dan de rest dunner gemaakt
for j=1:10
	AVER=thinnerwall(AVER,AITRI,endo,0.92,12,2.1,0);
	AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry	
	orgAVER=AVER;	
end
 
orgAVER=AVER;
disp('smooth heart');
AVER=smoothSurfaceVer(AVER,AITRI,endo);
AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry
orgAVER=AVER;
savetri([outname 'a.tri'],AVER/1000,AITRI);
TRIstatistics(AVER,AITRI,2000);


% Er zijn wat punten bliven hangen, die worden naar het gemiddelde van zijn
% buren getrokken. Er zitten ook wat epi vertices bij.
orgAVER=AVER;
nds=[34 788 796 908 1126 516 437 948 437 948 26 28 26 28 27 28	];
for i=1:length(nds)
	[a,b]=find(AITRI==nds(i)); 
	neigh=[unique(AITRI(a,:))']; neigh=neigh(neigh~=nds(i));
	AVER(nds(i),:)=mean(AVER(neigh,:));
end
AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry
orgAVER=AVER;

% nog een laatste slag
AVER=thinnerwall(AVER,AITRI,endo,0.86,12,2.1,0);
AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry
orgAVER=AVER;
AVER=thinnerwall(AVER,AITRI,endo,0.86,12,2.1,3);
AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry
orgAVER=AVER;

% orgAVER=AVER;

% disp('smooth heart');
% AVER=smoothSurfaceVer(AVER,AITRI,nds);
% AVER=CheckGeom(AVER,AITRI,orgAVER);%check geometry
% orgAVER=AVER;


TRIstatistics(AVER,AITRI,4000);

savetri([outname '.tri'],AVER/1000,AITRI);





% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VER=thinnerwall(VER,ITRI,nodes,delta,maxdist,mindist,MIND)

% VER vertices
% ITRI
% noded : the nodes (indices) that will be moved
% maxdist: The maximal wall distance for a node to be moved. Only nodes
% with a smaller wall thickness than max dist will be moved
% mindist: minamal wall thickness
% MIND: oposite to maxdist. The wall thickness must be greater than this
% distance to be changed. Thsi was needed for the left auricle

c=zeros(length(nodes),3);
for i=1:length(nodes)
	[ti,l]=find(ITRI==nodes(i));
	b=cross(VER(ITRI(ti(:),2),:)-VER(ITRI(ti(:),1),:),VER(ITRI(ti(:),3),:)-VER(ITRI(ti(:),1),:));
	c(i,:)=mean(b); c(i,:)=c(i,:)./norm(c(i,:)); % mm
end


MAXD=50;
k=0;kk=0;
for i=1:length(nodes)
	v2=VER(nodes(i),:)+c(i,:);
	TR=linetris(VER,ITRI,VER(nodes(i),:),v2);
	D=min(TR(TR(:,5)>0.1,5));
	if ~isempty(D) & D < maxdist & D >MIND
		keepVert=VER(nodes(i),:);
		alpha=max(mindist,D*delta);
		VER(nodes(i),:)=VER(nodes(i),:)+(D-alpha)*c(i,:);
		
		[ti,l]=find(ITRI==nodes(i));
		tj=unique(ITRI(ti,:));
		tj(find(tj==nodes(i)))=[];	k=k+1;
		for j=1:length(tj)
			TRT=linetrisect(VER,ITRI,VER(nodes(i),:),VER(tj(j),:));
			if ~isempty(TRT)
				keepVert;
				VER(nodes(i),:)=keepVert;  %restore due to crossection
				k=k-1;
				kk=kk+1;
				break;
			end
		end
	end
end

disp(['adapted ' num2str(k) ' nodes. Changed wall thickness  with a factor ' num2str(delta) ' maximal wall thickness ' num2str(maxdist) ' mm']);


%------------------------------------------------------------------

function otherVER=UpdateVertices(otherVER,VER,orgVER)
% adapt also the endocadrium, e.g. the vertices must be adapted

for i=1:length(VER)
	a=find(otherVER(:,1)==orgVER(i,1) & otherVER(:,2)==orgVER(i,2) & otherVER(:,3)==orgVER(i,3));
	if ~isempty(a)
		otherVER(a,:)=VER(i,:);
	end
end

%**************************************************************************
function VER=smoothSurfaceVer(VER,ITRI,vers)

hw=waitbar(0,'smoothing the surface');
for i=1:length(vers);
	keep=[1:length(VER)];
	keep(vers(i))=[]; 
	T=intripol(VER,ITRI,keep);
	VER=T*VER(keep,:);
	waitbar(i/length(vers),hw);
end
close(hw);


% %------------------------------------------------------------------
% 
% function plotAuto(AVER,AITRI)
% 
% 
% disp('plot auto solid angle')
% B=zeros(length(AVER));
% hw=waitbar(0,'Calculating auto solid angle');
% for i=1:length(AVER)
% 	[B(i,:),jsing]=rowforw(AVER,AITRI,AVER(i,:));
% 	waitbar(i/length(AVER),hw);
% end
% close(hw);
% AMA_A=B;
% 
% colA=AMA_A.*eye(size(AMA_A)); colA=colA(colA~=0);
% patch ('Faces',AITRI,'Vertices',AVER,'FaceVertexCData',colA,'FaceColor','interp',...
% 			'facealpha',1,'FaceLighting','phong','edgecolor','k','ButtonDown','selectnode');
% axis off normal;
% view(90,0); colormap(loadasci('tims.mcm'));
% %caxis([-max(abs(colA)) max(abs(colA))]);
% colorbar;axis equal vis3d
