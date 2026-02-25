% function [ITRI,ITRI_2,ITRI_3,ITRI_4,ITRI_5]=retriangulate(VER,ITRI,HVER,ITRI_1,VER_2,ITRI_2,VER_3,ITRI_3,VER_4,ITRI_4,VER_5,ITRI_5)
% scripts called: (included)
%                 [QUAD,tri1,tri2]=findQuad(VER,ITRI,vers1,vers2);
%                 [quadrang]=quadrangles_pvd(VER,edge)

function varargout=retriangulate(varargin)

if mod(length(varargin),2) ~= 0
	disp('Inputs should be paired: vertices and index (VER ITRI)');
	return;
elseif nargout ~=nargin/2
	disp('more inputs than outputs, only the indices are returned');
	return
else
	VER=varargin{1};
	ITRI=varargin{2};
	varargout{1}=ITRI;
	ngeoms=length(varargin)/2-1;
	for i=2:length(varargin)/2
		eval(['VER_' num2str(i-1) '= varargin{i*2-1};']);
		eval(['ITRI_' num2str(i-1) '= varargin{i*2};']);
		varargout{i}=eval(['ITRI_' num2str(i-1) ';']);
	end
end

% change triangles
disp('change triangles');
changed=100; nn=0;
totCh=0;

figure(33);clf; 
patch ('Faces',ITRI,'Vertices',VER,'FaceColor','r','Facealpha',1);
axis off equal vis3d; view(90,0)
hold on

while nn<10 & changed >2
	nn=nn+1;
	changed=0;
	hw=waitbar(0,['reorder triangulization, nr: ' num2str(nn)]);
	% determine all edges
	adj=tri2graph(VER,ITRI,0);
	adj=triu(adj);
	[a,b]=find(adj>0);
	edges=[a b];

   for ei=1:length(edges),
		vers1=VER(edges(ei,1),:);
		vers2=VER(edges(ei,2),:);
		[QUAD,t1,t2]=findQuad(VER,ITRI,vers1,vers2);
		
		if numel(QUAD)==5
			% try to improve the QUAD that consist of triangle t1 and t2
			angs=quadrangles_pvd(VER,QUAD);
			[ma ima]=max(angs);
			if (QUAD(ima)~=edges(ei,1) & QUAD(ima)~=edges(ei,2)) 
				% determine area ratio's of new and 'old' traingles
				a1=abs(det([ 1 1 1;VER(ITRI(t1,2),:)-VER(ITRI(t1,1),:);VER(ITRI(t1,3),:)-VER(ITRI(t1,1),:)]))./2;
				a2=abs(det([ 1 1 1;VER(ITRI(t2,2),:)-VER(ITRI(t2,1),:);VER(ITRI(t2,3),:)-VER(ITRI(t2,1),:)]))./2;
				deltO=abs((a1-a2)/(a1+a2));
				a1=abs(det([1 1 1;VER(QUAD(2),:)- VER(QUAD(1),:);VER(QUAD(4),:)-VER(QUAD(1),:)]))/2;
				a2=abs(det([1 1 1;VER(QUAD(2),:)- VER(QUAD(3),:);VER(QUAD(4),:)-VER(QUAD(3),:)]))/2;
				deltN=abs((a1-a2)/(a1+a2));
						  
				% determine angles beteen the new and old triangles
				n1=cross(VER(QUAD(2),:)-VER(QUAD(1),:),VER(QUAD(4),:)-VER(QUAD(1),:));n1=n1/norm3d(n1);
				n2=cross(VER(QUAD(2),:)-VER(QUAD(3),:),VER(QUAD(4),:)-VER(QUAD(3),:));n2=n2/norm3d(n2);
				nang=acos(sum(n1.*n2));
				nold1=cross(VER(QUAD(1),:)-VER(QUAD(2),:),VER(QUAD(3),:)-VER(QUAD(2),:));nold1=nold1/norm3d(nold1);
				nold2=cross(VER(QUAD(1),:)-VER(QUAD(4),:),VER(QUAD(3),:)-VER(QUAD(4),:));nold2=nold2/norm3d(nold2);
				nOldang=acos(sum(nold1.*nold2));
				% determin if the 'new' edge will crossect a triangle
				inter=linetrisect(VER,ITRI,VER(QUAD(2),:),VER(QUAD(4),:));
				
				% when no crossesction is found and the angles between new and
				% old  triangle sdiffer less than 30  and the area ratios
				% differ not more than 30% .... go on
				if isempty(inter(inter(:,5)>0 & inter(:,5)<1,:)) & ...
					 abs(nang-nOldang) < pi/6 & (deltN > deltO)
					% find here in the other geoms. If only 3 vertices are found in another
					% geometry, this quad cannot be changed (boundary)
					
					v1=QUAD(2);
					v2=QUAD(4);
					tri1=find((ITRI(:,1)==v1 & ITRI(:,2)==v2) |...
								 (ITRI(:,2)==v1 & ITRI(:,3)==v2) |...
								 (ITRI(:,3)==v1 & ITRI(:,1)==v2));
					tri2=find((ITRI(:,1)==v2 & ITRI(:,2)==v1) |...
								 (ITRI(:,2)==v2 & ITRI(:,3)==v1) |...
								 (ITRI(:,3)==v2 & ITRI(:,1)==v1));
					restore =0;
					if ~isempty(tri1) |~isempty(tri2)
						restore=1
						vers1=VER(v1,:);
						vers2=VER(v2,:);
						% restore here
						[QUAD,t1,t2]=findQuad(VER,ITRI,vers1,vers2);
					end
					% check for border or if thhis geometry intersects with any
					% of the other geometries
					border=0;
					intersection=0;
					for ng=1:ngeoms
						id=num2str(ng);
						inter=eval(['linetrisect(VER_' id ',ITRI_' id ',VER(QUAD(2),:),VER(QUAD(4),:));']);
						intersection=~isempty(inter(inter(:,5)>0 & inter(:,5)<1,:));
						if intersection
							break;
						end;
						eval(['[QUAD_' id ',t' id '_1,t' id '_2]=findQuad(VER_' id ',ITRI_' id ',vers1,vers2);']);
						if eval([' ~isempty(QUAD_' id ') '])
							if eval(['numel(QUAD_' id ')~=5  '])
								border=1;
								break;
							else
								for qi=1:length(QUAD)
									if eval(['isempty(find(VER_' id '(:,1)==VER(QUAD(qi),1) & VER_' id '(:,2)==VER(QUAD(qi),2) & VER_' id '(:,3)==VER(QUAD(qi),3)))'])
										border=1;
										break;
									end
								end
							end
						end
					end
					% if one of the other geometries contain only one tringle, it
					% does mean this edge is shared by more than one geometry. It
					% cannot be changed therefore.
					if ~border & ~intersection
						%'adapt'
						if (~restore) 
							changed=changed+1;
						end
						ITRI(t1,1:3)=[QUAD(1) QUAD(2) QUAD(4)];
						ITRI(t2,1:3)=[QUAD(3) QUAD(4) QUAD(2)];
						% draw the line
						line(VER(QUAD([2,4]),1), VER(QUAD([2,4]),2), VER(QUAD([2,4]),3),'col','w')
						for ng=1:ngeoms
							id=num2str(ng);
							if eval(['numel(QUAD_' id ')==5'])
								angs=eval(['quadrangles_pvd(VER_' id ',QUAD_' id '); ']);
								[ma ima]=max(angs);
								eval(['ITRI_' id '(t' id '_1,1:3)=[QUAD_' id '(1) QUAD_' id '(2) QUAD_' id '(4)];']);
								eval(['ITRI_' id '(t' id '_2,1:3)=[QUAD_' id '(3) QUAD_' id '(4) QUAD_' id '(2)];']);
							end
						end
					end
				end
			end
		end
		waitbar(ei/length(edges));
	end
	close(hw);
	totCh=totCh+changed;
	disp(['changed this loop: ' num2str(changed) '   round: ' num2str(nn) '  tot changed: ' num2str(totCh)]);
end

varargout{1}=ITRI;
for i=2:length(varargin)/2
	varargout{i}=eval(['ITRI_' num2str(i-1) ';']);
end


%---------------------------------------------------------------
function [QUAD,tri1,tri2]=findQuad(VER,ITRI,vers1,vers2);
v1=[];v2=[];j1=0;j2=0;
QUAD=[];tri1=0;tri2=0;

% vers1 and vers2 are one rib belonging to two triangles. Determine the
% QUAD of those two triangles
v1=find(VER(:,1)==vers1(1) & VER(:,2)==vers1(2) & VER(:,3)==vers1(3));
v2=find(VER(:,1)==vers2(1) & VER(:,2)==vers2(2) & VER(:,3)==vers2(3));
if isempty(v1) |isempty(v2)
	return;
end

tri1=find((ITRI(:,1)==v1 & ITRI(:,2)==v2) |...
			 (ITRI(:,2)==v1 & ITRI(:,3)==v2) |...
			 (ITRI(:,3)==v1 & ITRI(:,1)==v2));
tri2=find((ITRI(:,1)==v2 & ITRI(:,2)==v1) |...
			 (ITRI(:,2)==v2 & ITRI(:,3)==v1) |...
			 (ITRI(:,3)==v2 & ITRI(:,1)==v1));
if isempty(tri1) |isempty(tri2) 
	QUAD=1; % generates an error
	return;
end
if length(tri1)>1 | length(tri2)>1
	QUAD=1;
	tri1
	tri2
    return;
end

if tri1==tri2
	QUAD=1;
	tri1
	tri2
    return;
end


QUAD=[v1,...
		ITRI(tri2,find(ITRI(tri2,:)~=v1 & ITRI(tri2,:)~=v2)),...
		v2,...
		ITRI(tri1,find(ITRI(tri1,:)~=v1 & ITRI(tri1,:)~=v2)),...
		v1];
		 
		 
% quadrangles.m
% function [quadrang,edge]=(VER,ITRI)
% compute the (internal) angles at the four vertices of the 
% quadrangle formed by two connected triangles 
% ITRI specifies three vertices for each if the TWO triangles
% VER specifies the vertex coordinates
% see also triangles
% 050201 a van oosterom

function [quadrang]=quadrangles_pvd(VER,edge)
    R(1,:)=VER(edge(2),:)-VER(edge(1),:);
    R(2,:)=VER(edge(3),:)-VER(edge(2),:);
    R(3,:)=VER(edge(4),:)-VER(edge(3),:);
    R(4,:)=VER(edge(1),:)-VER(edge(4),:);
    lr=norm3d(R);
    for i=1:4,
        j=icyc(i-1,4);, a=R(i,:); b=-R(j,:);
        sharp=sign(a*b');
        la=lr(i); lb=lr(j);
        C(i,:)=cross(a,b);
        acrb=norm3d(cross(a,b));
        sinang=asin(max(min(acrb/(la*lb),1),-1));
        if sharp<0, sinang=pi-sinang; end
        quadrang(i)=sinang;
             
    end
  % treat possibility of angles > pi 
  scc=sum(sign(C*C'));
  [mi imi]=min(scc);
  if mi<0, quadrang(imi)=2*pi-quadrang(imi);end 