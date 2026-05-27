function varargout=calcAdjT(varargin) %[adjact,velocity,minvelocity,ambiPaths]
%% calulate the velocity on the heartsurface based on the given activation and adjacency matrix
%
% Create an adjacency matrix with the given depolarization moments as the 
% 'distances' with the restriction that the condution velocity to get from
% one node to the other is less than conduction velocity (optinal parameter 
% 3). 
% H.ADJ = 3D distances adjacency
% H.ADJW = ADJ2 as you just compute them
% H.DISTW = distance matrix base on H.ADJW

% output is a new scaled adjacency matrix that you ned to use. 


if length(varargin) < 2
	error('This routine need the depolarization moments (1) and the adajcency matrix (2)');
else
	dep=varargin{1};
	H=varargin{2};
	locMinArea=20;
	conduction_velocity=1; % m/s
	pp=3;
	while pp<=length(varargin)
		if ischar(varargin{pp})
			key=lower(varargin{pp});
			switch key
				case 'velocity'
					conduction_velocity=varargin{pp+1};pp=pp+2;
				otherwise
					error('unknown parameter');
			end			
		end
	end
	
end

%% Determine nodes at which the sinus wavefront (depV) has a local minimum (locmin)
locmin=[];
for i=1:length(dep)
	if all(dep(H.ADJ(:,i)>0& H.ADJ(:,i)<locMinArea)-dep(i)>=0)
		locmin=[locmin i];
	end
end

% Through the wall distances are twice as long, because the conduction
% velocity is slower transverse

Vdep=H.ADJ;
Tdep=zeros(size(Vdep));
for i=1:length(Vdep)
	dt=(dep-dep(i))';
	Vdep(i,dt>0)=Vdep(i,dt>0)./dt(dt>0);
	Tdep(i,:)=dt;
end
Vgeom=conduction_velocity*H.ADJ./H.ADJW;Vgeom(isnan(Vgeom))=0;
adjT=H.ADJW./conduction_velocity;
adjT(Vgeom>Vdep)=abs(Tdep(Vgeom>Vdep));


%% adapt adjT and determine the slowest velocities
deps=[(1:length(dep))' dep];	deps=sortrows(deps,2);
fixed=zeros(length(dep),1);		fixed(locmin)=1;
minvelocity=zeros(size(fixed));	minvelocity(locmin)=0; velocity=minvelocity;
kkk=0;

for i=1:length(fixed)
	if ~fixed(deps(i,1)) 
		af=find(fixed==1 & H.ADJ(:,deps(i,1))>0);
		% use the corrected adjacency values to mimic slower connections
		% through the wall
		mydist=H.DISTW(af,deps(i,1)); 
        mytime=deps(i,2)-dep(af);
        myvelo=mydist./mytime;
		% some nodes are fixed because they are foci. In some cases this
		% can lead to negative intervals
		af(myvelo<=0 | abs(myvelo)==Inf)=[]; 
		myvelo(myvelo<=0 | abs(myvelo)==Inf)=[];
		con=af(myvelo==min(myvelo)); 
		adjT(deps(i,1),con(1))=deps(i,2)-dep(con(1));
		adjT(con(1),deps(i,1))=adjT(deps(i,1),con(1));
		minvelocity(deps(i,1))=min(myvelo);
% 		geomVelo=Vgeom(deps(i,1),af);
% 		af=af(myvelo>0 & myvelo < geomVelo);				
% 		if min(myvelo) < conduction_velocity && ~isempty(af)
% 			adjT(deps(i,1),af)=abs(depV(af)-deps(i,2));
% 			adjT(af,deps(i,1))=adjT(deps(i,1),af);
% 			kkk=kkk+2*(length(af)-1);
% 		end
		fixed(deps(i,1))=1;
	end	
end

varargout{1}=adjT;
if nargout==1,	return;end
varargout{2}=locmin;
if nargout==2,	return; end
varargout{3}=minvelocity;
if nargout==3,	return;end




return
%% Check for the paths to determine the velocity

depsi=[];
for i=1:length(locmin)
	eval(['[depsi(i,:),path' num2str(i) ']=graphdistone(adjT,locmin(i));']);
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


