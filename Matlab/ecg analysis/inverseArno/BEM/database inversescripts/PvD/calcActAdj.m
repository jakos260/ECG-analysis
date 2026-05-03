function varargout=calcActAdj(varargin) %[adjact,velocity,minvelocity,ambiPaths]
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
	myover=varargin{3}; % vertices are clamped at the conduction_velocity value
	if isempty(myover)
		myover=zeros(size(depV));
	end
	conduction_velocity=1; % m/s
	myovelocity=conduction_velocity;
	purkinjeVelocity=3;	
	if length(varargin) > 2 
		if length(varargin) >= 4 && length(varargin{4})==1
			conduction_velocity=varargin{4};
			if length(varargin) >= 5 && length(varargin{5})==1
				myovelocity=varargin{5};
				if length(varargin) > 5
					VER =varargin{6};
					ITRI=varargin{7};
					if length(varargin) > 7
						disp('Parameters 7 and further are ignored')
					end
				end
			end
		elseif length(varargin) > 4
			VER =varargin{4};
			ITRI=varargin{5};
			if length(varargin) >4
				disp('Parameters 5 and further are ignored')
			end
		else 
			disp('Parameters 3 and further are ignored')
		end
	end
end
adj=VENTR.ADJ2W;
distadj=VENTR.ADJ;

%% Determine nodes at which the sinus wavefront (depV) has a local minimum (locmin)
locmin=calcLocMin();

dist=VENTR.ADJ;

Vdep=VENTR.ADJ;
Tdep=zeros(size(Vdep));
for i=1:length(Vdep)
	dt=abs(VENTR.depV-VENTR.depV(i))';
	Vdep(i,dt>0)=Vdep(i,dt>0)./dt(dt>0);
	Tdep(i,:)=dt;
end

Vpurkgeom=purkinjeVelocity*VENTR.ADJ./adj;Vpurkgeom(isnan(Vpurkgeom))=0;
adjact=adj./purkinjeVelocity;
adjact(Vpurkgeom>Vdep)=Tdep(Vpurkgeom>Vdep);
if ~isempty(find(myover==1))
	Vgeom=myovelocity*VENTR.ADJ./adj;Vgeom(isnan(Vgeom))=0;
	adjmyo=adj./myovelocity;
	adjmyo(Vgeom>Vdep)=abs(Tdep(Vgeom>Vdep));
	
	adjact(myover==1,:)=adjmyo(myover==1,:);
	adjact(:,myover==1)=adjmyo(:,myover==1);
end
varargout{1}=adjact;
if nargout==1
	return;
end
if nargout==2
	MV=VENTR.ADJ./adjact;MV(isnan(MV))=Inf;
	varargout{2}=min(MV);
	return;
end


%% adapt adjact and determine the slowest velocities
deps=[(1:length(depV))' depV];	deps=sortrows(deps,2);
fixed=zeros(size(adj,1),1);		fixed(locmin)=1;
minvelocity=zeros(size(fixed));	minvelocity(locmin)=0; velocity=minvelocity;
kkk=0;
for i=1:length(fixed)
	if ~fixed(deps(i,1)) 
		af=find(fixed==1 & adj(:,deps(i,1))>0);
		if isempty(af)
			continue;
		end;		
		% use the corrected adjjacency values to mimic slower connections
		% through the wall
		mydist=distadj(af,deps(i,1)); mytime=deps(i,2)-depV(af);myvelo=abs(mydist./mytime);
		% some nodes are fixed because they are foci. In some cases this
		% can lead to negative intervals
		if isempty(find(myvelo>0 & myvelo<=purkinjeVelocity))
			af(myvelo>min(myvelo))=[]; 
			myvelo(myvelo>min(myvelo))=[]; 
		else
			af(myvelo<=0 | myvelo>purkinjeVelocity)=[]; 
			myvelo(myvelo<=0 | abs(myvelo)>purkinjeVelocity)=[];
		end
		con=af(myvelo==min(myvelo)); 
		if isempty(con) || (myover(deps(i,1)) || myover(con(1)))
			continue;
		end;	
% 		adjact(deps(i,1),con(1))=deps(i,2)-depV(con(1));
% 		adjact(con(1),deps(i,1))=adjact(deps(i,1),con(1));
		minvelocity(deps(i,1))=min(myvelo);
% 		geomVelo=Vgeom(deps(i,1),af);
% 		af=af(myvelo>0 & myvelo < purkinjeVelocity);				
		if ~isempty(af)&& min(myvelo) < purkinjeVelocity 
			adjact(deps(i,1),af)=abs(depV(af)-deps(i,2));
			adjact(af,deps(i,1))=adjact(deps(i,1),af);
			kkk=kkk+2*(length(af)-1);
		else
			
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


