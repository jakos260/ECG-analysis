function GEOM=initfi(GEOM)
% index of the vertex of the right apex


% to make sure the contours around the valves are drawn properly
GEOM.endoVER(440)=0;    
GEOM.endoVER(462)=0;        
GEOM.endoVER(112)=0;            
GEOM.endoVER(439)=0;  
GEOM.endoVER(438)=0;             

GEOM.frankLeads=[235   322   177   165   272    82   162];
% GEOM.first6Finlay=[286    17   118   155   391   209];
GEOM.testLeads=loadmat('fi.leads');
GEOM.testLeadsrd=loadmat('fird.leads');
GEOM.v12=[12 18 25 32 40 45 63 64 65 ];
buur=graphdist(GEOM.ITRI);
D=GEOM.ADJ;D(buur==1)=0;
maxd=max(max(D)); % approximatly the distance between apex and base

for i=1:length(GEOM.endoVER)
	if any(GEOM.VER(i,1)==GEOM.RVER(:,1) & ...
		   GEOM.VER(i,2)==GEOM.RVER(:,2) & ...
		   GEOM.VER(i,3)==GEOM.RVER(:,3),1)
		GEOM.Lpurkinjever(i)=0;
	else
		GEOM.Rpurkinjever(i)=0;
	end
	if GEOM.Lpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd-35)*0.3
		GEOM.Lpurkinjever(i)=0;
	end
	if GEOM.Rpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd-25)*0.3
		GEOM.Rpurkinjever(i)=0;
	end
end