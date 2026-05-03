function varargout=calcActAdjLimited(varargin) %[adjact,velocity,minvelocity,ambiPaths]
%% calulate the velocity on the heartsurface based on the given activation and adjacency matrix
%
% Create an adjacency matrix with the given depolarization moments as the 
% 'distances' with the restriction that the condution velocity to get from
% one node to the other is less than conduction velocity (optinal parameter 
% 3). 

global VENTR
% global ATRIA

if length(varargin) < 2
	error('This routine need the depolarization moments (1) and the adajcency matrix (2)');
else
	depV=varargin{1};
	adj =varargin{2};
	locMinArea=20;
	conduction_velocity=1; % m/s
	if length(varargin) > 2 
		if length(varargin{3})==1
			conduction_velocity=varargin{3};
			if length(varargin) > 4
				VER =varargin{4};
				ITRI=varargin{5};
				if length(varargin) > 5
					disp('Parameters 6 and further are ignored')
				end
			end
		elseif length(varargin) > 3
			VER =varargin{3};
			ITRI=varargin{4};
			if length(varargin) >4
				disp('Parameters 5 and further are ignored')
			end
		else 
			disp('Parameters 3 and further are ignored')
		end
	end
end

distadj=adj;

%% Determine nodes at which the sinus wavefront (depV) has a local minimum (locmin)
locmin=[];
for i=1:length(depV)
	if all(depV(adj(:,i)~=0 & adj(:,i)<locMinArea)-depV(i)>0)
		locmin=[locmin i];
	end
end

dist=adj;
% 	Through the wall distances are twice as long, because the conduction velocity is slower transverse
% if exist('VENTR','var') && isfield(VENTR,'ITRI') 
% 	adj(~buur)=adj(~buur)*2;
% end
if exist('VENTR','var') && isfield(VENTR,'ADJW')
	adj=VENTR.ADJW;
elseif exist('VENTR','var') && isfield(VENTR,'ITRI')
	buur=graphdist(VENTR.ITRI);
	adj(~buur)=adj(~buur)*2;
elseif exist('ITRI','var') 
	buur=graphdist(ITRI);
	adj(~buur)=adj(~buur)*2;	
end
% if exist('ATRIA','var') && isfield(ATRIA,'endoVER')
% 	for i=1:length(adj) 
% 		adj(ATRIA.endoVER~=ATRIA.endoVER(i),i)=adj(ATRIA.endoVER~=ATRIA.endoVER(i),i)*2;
% 		adj(i,ATRIA.endoVER~=ATRIA.endoVER(i))=adj(i,ATRIA.endoVER~=ATRIA.endoVER(i))*2;		
% 	end
% end


adjact=adj./conduction_velocity;

Vdep=VENTR.ADJ;
Tdep=zeros(size(Vdep));
for i=1:length(Vdep)
	dt=(VENTR.depV-VENTR.depV(i))';
	Vdep(i,dt>0)=Vdep(i,dt>0)./dt(dt>0);
	Tdep(i,:)=dt;
end
Vgeom=conduction_velocity*VENTR.ADJ./VENTR.ADJW;Vgeom(isnan(Vgeom))=0;
adjact=VENTR.ADJW./maxvelo;
adjact(Vgeom>Vdep)=abs(Tdep(Vgeom>Vdep));




%% adapt adjact and determine the slowest velocities
deps=[(1:length(depV))' depV];	deps=sortrows(deps,2);
fixed=zeros(size(adj,1),1);		fixed(locmin)=1;
minvelocity=zeros(size(fixed));	minvelocity(locmin)=0; velocity=minvelocity;
kkk=0;
for i=1:length(fixed)
	if ~fixed(deps(i,1)) 
		if deps(i,1)==1
			stop=1;
		end
		af=find(fixed==1 & adj(:,deps(i,1))>0);
		% use teh corrected adjjacency values to mimic slower connections
		% through the wall
		mydist=distadj(af,deps(i,1)); mytime=deps(i,2)-depV(af);myvelo=mydist./mytime;
		% some nodes are fixed because they are foci. In some cases this
		% can lead to negative intervals
		af(myvelo<=0 | abs(myvelo)==Inf)=[]; 
		myvelo(myvelo<=0 | abs(myvelo)==Inf)=[];
		con=af(myvelo==min(myvelo)); 
		if min(myvelo)<=conduction_velocity
			adjact(deps(i,1),con(1))=deps(i,2)-depV(con(1));
			adjact(con(1),deps(i,1))=adjact(deps(i,1),con(1));
		end
		minvelocity(deps(i,1))=min(myvelo);
% 		if minvelocity(deps(i,1)) > conduction_vel
% 		endocity
% 			stop=1;
		af=af(myvelo>0 & myvelo < conduction_velocity);				
		if min(myvelo) < conduction_velocity && ~isempty(af)
			adjact(deps(i,1),af)=abs(depV(af)-deps(i,2));
			adjact(af,deps(i,1))=adjact(deps(i,1),af);
			kkk=kkk+2*(length(af)-1);
		end
		fixed(deps(i,1))=1;
	end	
end






varargout{1}=adjact;
if nargout==1
	return;
end
if nargout==2
	varargout{2}=minvelocity;
	return;
end

%% Check for the paths to determine the velocity

depsi=[];
for i=1:length(locmin)
	eval(['[depsi(i,:),path' num2str(i) ']=graphdistone(adjact,locmin(i));']);
	depsi(i,:)=depsi(i,:)+depV(locmin(i));
end
if length(locmin)>1
	depsin=min(depsi)';
else
	depsin=depsi';
end
ambiPaths=zeros(size(depsin));
for i=1:length(depsin)
	if isempty(find(i==locmin,1))
		a=find(abs(depsi(:,i)-depsin(i))<1e-12);

		other=[];
		for j=1:length(a)
			eval(['other(j)=path' num2str(a(j)) '(i,1);']);
		end
		if length(a) >1 && ~isempty(find(other(1)~=other,1))
			ambiPaths(i)=1;%ambiPaths+1;
		end
		velocity(i)=min(dist(other,i)./abs(depV(i)-depV(other)));	
	end
end
if nargout>1
	varargout{2}=velocity;
end
if nargout>2
	varargout{3}=minvelocity;
end
if nargout>3
	varargout{4}=depsin;
end

if nargout>4
	varargout{5}=ambiPaths;
end

%% messages
disp(['Reduced the velocity for ' num2str(kkk/2+1500-length(locmin))  ' connections (='	...
	   num2str(100*kkk/length(nonzeros(adj)),3) ' %)   Default conduction velocity ' num2str(conduction_velocity) ' m/s']);
disp(['The mean propagation velocity ' num2str(mean(velocity),3) ' +/- ' num2str(std(velocity),3) ' m/s']);
if exist('ambiPaths','var')
	disp(['Number of abiguous routes: ' num2str(sum(ambiPaths))]); 
end
notequal=find(abs(depsin-depV)>1e-12);
if isempty(notequal)
	disp('All reconstructed activation moments of the nodes are equal to the given values')
else
	disp([ num2str(length(notequal)) ' nodes differ'])
end

%% show results
if exist('VER','var')
	if ~isempty(notequal) 
		figure(1);	showAtria(VER,ITRI,depsin-depV)

		figure(2);clf;colormap(loadmat('tims.mcm'))
		patch('Faces',ITRI,'Vertices',VER,...
			  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
			  'FaceVertexCData',depsin,'FaceColor','interp',...
			  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
		axis equal off;view(90,0);colorbar
	end
	figure(3);	clf;	colormap(loadmat('tims.mcm'))
	patch('Faces',ITRI,'Vertices',VER,...
		  'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,...
		  'FaceVertexCData',velocity,'FaceColor','interp',...
		  'edgecolor','none','FaceAlpha',1,'buttondownFcn','selectnode');
	axis equal off;view(90,0);colorbar
end


