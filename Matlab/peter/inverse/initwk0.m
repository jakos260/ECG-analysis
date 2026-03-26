function GEOM=initwk0(GEOM)
if ~strfind(GEOM.type,'atria')
	% index of the vertex of the right apex
	GEOM.IRAPEX=39;

	% to make sure the contours around the valves are drawn properly
	GEOM.endoVER(383)=0;
	GEOM.endoVER(384)=0;
	GEOM.endoVER(540)=0;    

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


		if GEOM.Lpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd-25)*0.4
			GEOM.Lpurkinjever(i)=0;
		end
		if GEOM.Rpurkinjever(i) && min(GEOM.DISTsurf(i,GEOM.endoVER==0))<=(maxd)*0.4
			GEOM.Rpurkinjever(i)=0;
		end			
	end
end
GEOM.frankLeads=[136   137   220   266   250    66   281];
% GEOM.first6Finlay=[261   121    22   171   160    72];
GEOM.testLeads=loadmat('wk0.leads');
GEOM.testLeadsrd=loadmat('wk0rd.leads');
