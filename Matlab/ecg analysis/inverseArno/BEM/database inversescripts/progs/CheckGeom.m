function varargout=CheckGeom(varargin)

% This program checks for intersections in a triangulated surface. The 
% program is meant to function after some vertices have been moved. If
% crossections exist the vertex is moved in the direction of the normal.
% Assuming the normal points outward for a closed surface most
% crossections will be solved. Finally the triangulated surface is
% checked again. If still crossections exist the vertex is moved to its
% original position (orgVER). Additionally the crossection between the
% surface at hand and another surface can be checked, for instance the
% blood cavity when the heart is moved.
% 
% INPUTS: 
% VER, ITRI, 
% orgVER (the positions of the vertices before any moving% of the vertices.
% VER_1,ITRI_1......VER_n,ITRI_n (other triangulated geometries)
% 
% Peter van Dam
% 12 April 2005
% 
%**************************************************************************

if mod(length(varargin)-3,2) ~= 0
	disp('Input: VER, ITRI,orgVER, VER1,ITRI1, .....VERn,ITRIn');
	return;
elseif nargout > 2 
	disp('outputs: newVER and error status');
	return
else
	VER=varargin{1};
	ITRI=varargin{2};
	orgVER=varargin{3};
	ngeoms=(length(varargin)-3)/2;
	for i=2:length(varargin)/2
		eval(['VER_' num2str(i-1) '= varargin{i*2};']);
		eval(['ITRI_' num2str(i-1) '= varargin{i*2+1};']);
	end
end

adj=tri2graph(VER,ITRI,0);
adj=triu(adj);
[a,b]=find(adj>0);
edges=[a b];
ner=0;error=1;
while error>0 & ner<=5
	% make sure the endo and epi do not intersect. If they do move away in
	% the direction of the normal, which should be pointed outward
	error=0;	adapted=0;
	hw=waitbar(0,'adapting crossections');
	for edg=1:length(edges)
		inter=linetrisect(VER,ITRI,VER(edges(edg,1),:),VER(edges(edg,2),:));
		nfin=1;
		while ~isempty(inter) & nfin<3 % stop if no crossing is found
			[ti1,l]=find(ITRI==edges(edg,1));
			[ti2,l]=find(ITRI==edges(edg,2));
			c1=mean(cross(VER(ITRI(ti1(:),1),:)-VER(ITRI(ti1(:),2),:),...
							  VER(ITRI(ti1(:),1),:)-VER(ITRI(ti1(:),3),:)));
			c2=mean(cross(VER(ITRI(ti2(:),1),:)-VER(ITRI(ti2(:),2),:),...
							  VER(ITRI(ti2(:),1),:)-VER(ITRI(ti2(:),3),:)));						 
			VER(edges(edg,1),:)=VER(edges(edg,1),:)-(nfin.*c1./norm(c1)); % move nfin mm outside
			VER(edges(edg,2),:)=VER(edges(edg,2),:)-(nfin.*c2./norm(c2)); % move nfin mm outside\
			inter=linetrisect(VER,ITRI,VER(edges(edg,1),:),VER(edges(edg,2),:));
% 			disp(['Adapted nodes of triangle ', num2str(tri) '  ' num2str(ITRI(tri,icyc(edg,3))), ' and ' num2str(ITRI(tri,icyc(edg+1,3))) ]);		
			adapted=adapted+1;
			nfin=nfin+1;
			error=1;
		end
		waitbar(edg/length(edges),hw);
	end
	close(hw);
	disp(['number adapted vertices: ' num2str(adapted)]);
	intersect=1;
	if error==0,intersect=0;end
	error=0;
	restoreI=[];
	for i=1:ngeoms
		id=num2str(i);
		eval(['VER_' id 't=UpdateVertices(VER_' id ',VER,orgVER);']); % temporary
	end

	hw=waitbar(0,['verifying crossections ']);
	for edg=1:length(edges)
		if intersect
			inter=linetrisect(VER,ITRI,VER(edges(edg,1),:),VER(edges(edg,2),:));
			if ~isempty(inter)
				restoreI=[restoreI; edges(edg,:)'];
% 				disp(['restored nodes of triangle ', num2str(tri) '  ' num2str(ITRI(tri,icyc(edg,3))), ' and ' num2str(ITRI(tri,icyc(edg+1,3)))]);			
				error=1;
			end
		end
		for i=1:ngeoms
			id=num2str(i);
			eval(['inter=linetrisect(VER_' id 't,ITRI_' id ',VER(edges(edg,1),:),VER(edges(edg,2),:));']);
			if ~isempty(inter)
				eval(['tri=ITRI_' id '(inter(1),:);'])
				f=[];
				for kk=1:3
					k=tri(kk);
					eval(['f=[f find(VER_' id '(k,1)==orgVER(:,1) & VER_' id '(k,2)==orgVER(:,2) & VER_' id '(k,3)==orgVER(:,3))];'])
				end
				restoreI=[restoreI edges(edg,:) f];
% 				disp(['restored nodes of triangle (other crossection)', num2str(tri) '  ' num2str(ITRI(tri,icyc(edg,3))), ' and ' num2str(ITRI(tri,icyc(edg+1,3)))]);			
				error=1;
			end
		end
		waitbar(edg/length(edges),hw);
	end
	close(hw);
	restoreI=unique(sort(restoreI));
	VER(restoreI,:)=orgVER(restoreI,:);
	disp(['number restored vertices: ' num2str(size(restoreI,1))]);
	ner=ner+1;
end
if error
	disp(['fatal error!!!!  After ' num2str(ner) ' runs, still crossections found.'])
end

varargout{1}=VER;
if nargout==2
	varargout{2}=error;
end

%------------------------------------------------------------------


function otherVER=UpdateVertices(otherVER,VER,orgVER)
% adapt also the endocadrium, e.g. the vertices must be adapted

for i=1:length(VER)
	a=find(otherVER(:,1)==orgVER(i,1) & otherVER(:,2)==orgVER(i,2) & otherVER(:,3)==orgVER(i,3));
	if a==959
		stop=1;
	end
	if ~isempty(a)
		otherVER(a,:)=VER(i,:);
	end
end


% 		finter=inter(inter(:,5)>0 & inter(:,5)<1,:);