function varargout=MIVesselActivation(varargin) 

global VENTR
vesselAreas=loadasci('C:\ECG_simulation\documents\articles\shortestpath\matlab\purkinje\vessels.asc');
if length(varargin) < 4
	error('This routine need the depolarization moments (1) and the adajcency matrix (2)');
else
	vesselArea=varargin{1};
	depV=varargin{2};
	repV=varargin{3};
	cvFactor=0;
	APDshort=0;
	transmural=1;
	blockFoci=0;
	repMethod=1;
	velocity=1;%ms-1
	purkinjevelocity=2;%ms-1
	purkinjever=[];
	remfocus=[];
	pSV=VENTR.pSV;
	restpot=1;
% 	defaults='geometry';
	defaults='myocard';
	doHighVelocities=1;
	ADJW=VENTR.ADJ;
	pp=4;
	while pp<=nargin
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'pv'
					cvFactor=varargin{pp+1};pp=pp+2;
				case 'apds'
					APDshort=varargin{pp+1};pp=pp+2;
				case 'trans'
					transmural=varargin{pp+1};pp=pp+2;
				case 'blockfoci'
					blockFoci=varargin{pp+1};pp=pp+2;
				case 'remfocus'
					remfocus=varargin{pp+1};pp=pp+2;
				case 'repmethod'
					repMethod=varargin{pp+1};pp=pp+2;
				case 'velocitydefaults'
					defaults=varargin{pp+1};pp=pp+2;
				case 'velocity'
					velocity=varargin{pp+1};pp=pp+2;
				case 'purkinjevelocity'
					purkinjevelocity=varargin{pp+1};pp=pp+2;
				case 'purkinjever'
					purkinjever=varargin{pp+1};pp=pp+2;					
				case 'up'
					pSV(1)=varargin{pp+1};pp=pp+2;
				case 'down'
					pSV(2)=varargin{pp+1};pp=pp+2;
				case 'tdomup'
					pSV(3)=varargin{pp+1};pp=pp+2;
				case 'tdomdown'
					pSV(4)=varargin{pp+1};pp=pp+2;
				case 'offset'
					pSV(5)=varargin{pp+1};pp=pp+2;
				case 'rest'
					restpot=varargin{pp+1};pp=pp+2;		
				case 'vessels'
					vesselArea=varargin{pp+1};pp=pp+2;		
				case 'adj'
					ADJW=varargin{pp+1};pp=pp+2;
					doHighVelocities=0;
				otherwise
					error('unknown parameter');
			end				
		else
			error('unknown parameter');
		end
	end	
	locMinArea=20;
end
if length(APDshort)>1 && length(APDshort)~=length(node)
	warning('APD shortening settings differs from the number ischemic areas only the first one is used');
	APDshort=APDshort(1);
end
if length(transmural)>1 && length(transmural)~=length(node)
	warning('Number of the transmural parameters values differs from the number ischemic areas only the first one is used');
	transmural=transmural(1);
end
if length(cvFactor)>1 && length(cvFactor)~=length(node)
	warning('Number of the propgation velocity factor parameters values differs from the number ischemic areas only the first one is used');
	cvFactor=cvFactor(1);
end

if ~isfield(VENTR,'A2')
	[VENTR.A2,VENTR.D2]=graphdist(VENTR.ITRI,VENTR.VER,4);
end
%%
if blockFoci
	blocktext=[' purkinje velocity reset to ' defaults];
else
	blocktext=[' purkinje velocity kept'];
end
if transmural==0
	disp(['Area ' num2str(vesselArea) ' velocity reduction ' num2str(cvFactor) '  APD shortning ' num2str(APDshort) ' single card only' '   velocity reset to: ' defaults blocktext]);
else
	disp(['Area ' num2str(vesselArea) ' velocity reduction ' num2str(cvFactor) '  APD shortning ' num2str(APDshort) ' transmural (factor ' num2str(transmural) ')' '   velocity reset to: ' defaults blocktext]);
end


%%
cvF=cvFactor;
for i=1:length(cvFactor)
	if cvFactor(i)<=-1,	cvF(i)=-1+1e-5;else cvF(i)=cvFactor(i);end
end
%%
% Find the nodes that are in the "ischemic area
imi=[];
if ~isempty(vesselArea)
	if transmural~=1	
		for k=1:length(vesselArea)
			a=find(vesselAreas==vesselArea(k));
			imi=[imi;a(VENTR.endoVER(a)==abs(transmural)) ];
		end
	else
		for k=1:length(vesselArea)
			imi=[imi;find(vesselAreas==vesselArea(k))];
		end
	end
	imi=unique(imi);
end

%% Determine nodes at which the sinus wavefront (depV) has a local minimum (locmin)
% A right ventricular block occurs in case the LAD is proximal occluced
% (vesselAreas =3

locmin=calcLocMin;
%% remlocmin
rems=zeros(size(locmin));
for i=1:length(locmin)
	if any(remfocus==locmin(i))
		rems(i)=1;
	end
end
remloc=locmin(rems==1);
if blockFoci
	for i=1:length(locmin)
% 		if vesselAreas(locmin(i))>=5 & vesselAreas(locmin(i))<=7 & any(vesselArea==3)
% 			rems(i)=1;
% else
		if any(imi==locmin(i),1)
			rems(i)=1;
		end
	end
	remloc=locmin(rems==1);
end

%%

if doHighVelocities % default
	adjactmi=calcActAdj(VENTR.depV,VENTR.ADJW,purkinjever,velocity);
	if 0
	Vdep=VENTR.ADJ;
	Tdep=zeros(size(Vdep));
	for i=1:length(Vdep)
		dt=abs(VENTR.depV-VENTR.depV(i))';
		Vdep(i,dt>0)=Vdep(i,dt>0)./dt(dt>0);
		Tdep(i,:)=dt;
	end
	Vgeom=velocity*VENTR.ADJ./ADJW;Vgeom(isnan(Vgeom))=0;
	Vgeommyo=Vgeom;
	if isempty(purkinjever)
		purkinjever=VENTR.endoVER;
		for i=1:length(purkinjever)
			if purkinjever(i) & min(VENTR.D2(i,VENTR.endoVER==0))<=25
				purkinjever(i)=0;
			end
			if purkinjever(i) & min(VENTR.D2(i,VENTR.endoVER==0))<=45
				purkinjever(i)=0;
			end
		end
	end
	for i=1:length(VENTR.ADJ),
		for j=1:length(VENTR.ADJ),
			if VENTR.A2(i,j)>0 && purkinjever(i) && purkinjever(j)
				Vgeom(i,j)=purkinjevelocity;	
			end; 
		end; 
	end
	adjactmi=VENTR.ADJ./velocity;
	adjactmi(Vgeom>Vdep)=abs(Tdep(Vgeom>Vdep));
	if ~isfield(VENTR,'orgs')
	d=1000000*ones(length(locmin),length(depV));
	for i=1:length(locmin)
		d(i,:)=graphdistone(adjactmi,locmin(i));
		d(i,:)=d(i,:)+depV(locmin(i));
	end
	md=min(d);
	orgs=zeros(length(md),length(locmin));
	for i=1:length(d)
		a=find(i==locmin);
		if isempty(a)
			[a,b]=find(md(i)==d);
		end
		orgs(i,a)=1;	
	end	
	end
	for i=1:length(rems)
		if rems(i)==1
			a=find(orgs(:,i)==1);
			for j=1:length(a)
				b=find(Vgeom(:,a(j))==purkinjevelocity);
				Vgeom(b,a(j))=Vgeommyo(b,a(j));
				Vgeom(a(j),b)=Vgeom(b,a(j));
			end
		end
	end
	adjactmi=VENTR.ADJ./velocity;
	adjactmi(Vgeom>Vdep)=abs(Tdep(Vgeom>Vdep));
	end
else				% in case an adjecancy matrix is given
	Vdep=VENTR.ADJ;
	Tdep=zeros(size(Vdep));
	for i=1:length(Vdep)
		dt=abs(VENTR.depV-VENTR.depV(i))';
		Vdep(i,dt>0)=Vdep(i,dt>0)./dt(dt>0);
		Tdep(i,:)=dt;
	end
	Vgeom=velocity*ones(size(ADJW));	
% 	for i=1:length(Vgeom)
% 		if VENTR.endoVER(i)
% 	        Vgeom(i,VENTR.endoVER==1)= Vgeom(i,VENTR.endoVER==1)*1; % purkinje velocities over the wall 2 m/s
% 		end
% 	end
	adjactmi=ADJW./velocity;
	adjactmi(Vgeom>Vdep)=abs(Tdep(Vgeom>Vdep));
end



%% Create adjacency matrix for the ischemic area

if cvF ~=0
	% reset only the adjacency matrix to the the geomtrical given values
	% when the velocity is higher than 1 ms-1 and the foci have to be
	% blocked. Thus the purkinje system is deactivated.
	if strcmp(defaults,'geometry')
		adjv=ADJW; % reset the adjacency matrix to the geomtrical given values
		for i=1:length(adjactmi)
			adjactmi(i,imi)=adjv(i,imi)./velocity;
			adjactmi(imi,i)=adjv(imi,i)./velocity;
		end
	else
		if blockFoci
			V=VENTR.ADJ./adjactmi;
			V(isnan(V))=0;
			adjactmi(V(imi,:)>velocity)=ADJW(V(imi,:)>velocity)./velocity;
			adjactmi(V(:,imi)>velocity)=ADJW(V(:,imi)>velocity)./velocity;			
		end
	end
	if length(cvF)==1
		timi=imi;
% 		if ~blockFoci
% 			for i=1:length(purkinjever)
% 				timi(timi==purkinjever(i))=[];
% 			end
% 		end
		adjactmi(timi,:)=adjactmi(timi,:)./(1+cvF);
		adjactmi(:,timi)=adjactmi(:,timi)./(1+cvF);
	else
		for i=1:length(cvF)
			eval(['timi=imi_' num2str(i)]);
			adjactmi(imi,:)=adjactmi(timi,:)./(1+cvF(i));
			adjactmi(:,imi)=adjactmi(:,timi)./(1+cvF(i));	
		end
	end
end

%% Determine ischemic deplorarization 
% Determine the new depolarization moments The local minima moments are
% also delayed (factor cvF)in case these nodes are within the ischemic
% area. When the ischemica area blocks the timing of the area is set to
% maximum value in depp.
% if cvF ~=0
	depmi=1000000*ones(length(locmin),length(depV));
	for i=1:length(locmin)
		if isempty(find(remloc==locmin(i),1)) %|| ~blockFoci 
			depmi(i,:)=	graphdistone(adjactmi,locmin(i))+depV(locmin(i));
			disp([num2str([i locmin(i) ]) '  used' ])
		else
			disp([num2str([i locmin(i) ]) '  removed' ])
		end
	end
	if length(locmin)>1
		depp=min(depmi)';
	else
		depp=depp';
	end
% else
% 	depp=depV;
% 	disp(['Kept original depolarzation sequence' ])
% end
% Align the ECGs, i.e. all depolarizations starts at min(depV) only use
% when the ventricular complex is simulated without P wave
% depp=depp-abs(min(depV)-min(depp));
% if abs(min(depV)-min(depp))~=0
% 	disp(['allignment ' num2str(abs(min(depV)-min(depp)))])
% end

% ADP within the ischemic area (the imi nodes) is shortend by the ADPshort 
% interval. The repolarization momnts are now the depp + ADP
tv=0:550;
if ~isempty(imi)
	if cvFactor<=-1 % colmplete block
		if repMethod~=1
			depp(imi)=max(depV)*6.01;
			depp(imi)=max(depp)*1.3+depV(imi);
			repp=repV;
			APD=repp-depV;
			repp(imi)=depp(imi)+APD(imi);
			tv=0:550+ceil(max(depp));
		else
			depp(imi)=max(depV)*1.01;
			depp(imi)=max(depp)*1.01;
			repp=repV;
			repp(imi)=max(repV)*1.01;
			repp(imi)=max(repp)*1.01;
		end
	else
		APD=repV-depV;
		fact=calcFact(VENTR.DIST,imi);
		APD(imi)=APD(imi)-(APDshort(min(length(APDshort),i))*(fact));
		repp=depp+APD;
	end
else
	repp=repV;
end
varargout{1}=depp;
if nargout>1
	varargout{2}=repp;
end
if nargout>2
	varargout{3}=imi;
end

if nargout>3	
	Tv=ones(length(depp),1)*tv;
	% Scale S between -0.9 and 0.1
	Svmi=gets_v(Tv,depp,repp,pSV,3);
	Svmi=Svmi-0.9;			
	if cvFactor>-1
		fact=calcFact(VENTR.DIST,imi);
		Svmi(imi,:)=Svmi(imi,:).*(1-((1-restpot).*(fact*ones(1,size(Svmi,2)))));
	else% colmplete block
		Svmi(imi,:)=0.1;
	end
	varargout{4}=Svmi;
	if nargout>4
		varargout{5}=Tv;
	end
end


%%
% clf
% fact=1./(1+exp(0.6*(VENTR.DIST(node,imi)-MIarea+6)));
% fact1=1-(VENTR.DIST(node,imi)./MIarea).^3;
% plot(VENTR.DIST(node,imi),fact,'b.',VENTR.DIST(node,imi),fact1,'r.')

function fact=calcFact(dists,imi)

fact=ones(size(imi));
nonischemicnodes=ones(size(dists,1),1);
nonischemicnodes(imi)=[];
for i=1:length(imi)
	fact(i)=min(dists(imi(i),nonischemicnodes));	
end
% fact=sort(fact);
fact=1./(1+exp((10-fact)));
% figure(1)
% plot(fact)
% figure(2)
% plot(fact1,'r')
% fact=1./(1+exp(0.6*(dists-radius+6)));
% fact=1./(1+exp(0.6*(dists-radius+6)));
% fact1=1-(VENTR.DIST(node,imi)./MIarea).^3;


