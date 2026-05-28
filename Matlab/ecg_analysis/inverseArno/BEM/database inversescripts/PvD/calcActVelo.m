function velo=calcActVelo(depV,DIST,maxdist)


velo=zeros(size(DIST,1),1);
for i=1:size(DIST,1),    
    a= find(DIST(i,:)>0 & DIST(i,:)<maxdist);
    velos=DIST(a,i)./(depV(a)-depV(i));
    % 
%     minv =min(velos(velos>0)); 
    minv =min(abs(velos));
    
    if~isempty(minv),
        velo(i)=minv;
    end;
end
return


%% calulate the velocity on the heartsurface based on the given activation and adjacency matrix
%
% Create an adjacency matrix with the given depolarization moments as the 
% 'distances' with the restriction that the condution velocity to get from
% one node to the other is less than conduction velocity (optinal parameter 
% 3). 


if nargin<3
	error('This routine need the depolarization moments (1) and the adajcency matrix (2) and the distance matrix');
elseif nargin==3
    locmin = calcLocmin(DIST,depV);
end
distadj=adj;

%% Determine nodes at which the sinus wavefront (depV) has a local minimum (locmin)
%% adapt adjact and determine the slowest velocities
deps=[(1:length(depV))' depV];	deps=sortrows(deps,2);
fixed=zeros(size(adj,1),1);		fixed(locmin)=1;
minvelocity=zeros(size(fixed));	
minvelocity(locmin)=0; 
velocity=minvelocity;
for i=1:length(fixed)
	if ~fixed(deps(i,1)) 
		af=find(fixed==1 & adj(:,deps(i,1))>0);
		if isempty(af)
			continue;
		end;		
		% use the corrected adjjacency values to mimic slower connections
		% through the wall
		mydist=distadj(af,deps(i,1)); 
        mytime=deps(i,2)-depV(af);
        myvelo=mydist./mytime;
		% some nodes are fixed because they are foci. In some cases this
		% can lead to negative intervals
		af(myvelo<=0 | abs(myvelo)==Inf)=[]; 
		myvelo(myvelo<=0 | abs(myvelo)==Inf)=[];
		con=af(myvelo==min(myvelo)); 
		if isempty(con)
			continue;
		end;	
		minvelocity(deps(i,1))=min(myvelo);
		fixed(deps(i,1))=1;
	end	
end
for i=1:length(locmin)
    minvelocity(locmin(i))= max(minvelocity(adj(locmin(i),:) > 0));
end

return
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


